import subprocess
import csv
from pathlib import Path
import pandas as pd
import Levenshtein
from Bio import SeqIO
from Bio.Seq import Seq

probescrape_shorter_list = {}
with open('probescrape-results.parsed.csv') as probescrape_results_csv:
    for line_dict in csv.DictReader(probescrape_results_csv):
        shorter_list_key = f"{line_dict['probe_seq']}_{line_dict['fw_seq']}_{line_dict['rev_seq']}"
        if shorter_list_key not in probescrape_shorter_list: probescrape_shorter_list[shorter_list_key] = line_dict

probescrape_shorter_list: list = [result_dict for result_dict in probescrape_shorter_list.values()]
def _output_results(_results: list, _name_str: str) -> None: pd.DataFrame(_results).to_csv(f'{_name_str}', index=False); return None
def _check_specificity(_probe_seq: str, _forward_primer: str, _reverse_primer: str, _targets_list: list, _non_targets_list: list) -> tuple:
    """
    """
    with open('primer', 'w') as primer_file: primer_file.write(f'M_BOVIS\t{_forward_primer}\t{_reverse_primer}\t50\t1000\n')

    result = subprocess.run([
        'ipcress',
        '-i', 'primer',
        '--sequence', 'db_targets-and-non-targets.fna',
        '-m', '1',
        '-p', 'FALSE',
        '-P', 'TRUE',
        ], capture_output=True)
    with open('products.fasta', 'w') as products_file:
        for line in result.stdout.decode().split('\n'):
            if 'ipcress' in line: continue
            if line.startswith('>'): line = f">{line.strip().split('seq ')[1]}"
            products_file.write(line.strip()+'\n')

        target_probe_binds: int = 0
        non_target_probe_binds: int = 0
        target_primer_score: float = 0
        non_target_primer_score: float = 0
        for seq_object in SeqIO.parse('products.fasta', 'fasta'):
            f_lev_distance: int = Levenshtein.distance(_forward_primer, seq_object.seq[:len(_forward_primer)])
            r_lev_distance: int = Levenshtein.distance(Seq(_reverse_primer).reverse_complement(), seq_object.seq[-len(_forward_primer):])
            primer_score = (f_lev_distance+r_lev_distance)/(len(_forward_primer)+len(_reverse_primer))
            try: accession = '_'.join(seq_object.id.split('_')[:seq_object.id.split('_').index('GCA')+2])
            except: accession = '_'.join(seq_object.id.split('_')[:seq_object.id.split('_').index('GCF')+2])
            if accession in _targets_list: 
                if _probe_seq in seq_object.seq: target_probe_binds += 1
                target_primer_score += primer_score
            elif accession in _non_targets_list:
                if _probe_seq in seq_object.seq: non_target_probe_binds += 1
                non_target_primer_score += primer_score

    print(target_probe_binds, non_target_probe_binds, target_primer_score, non_target_primer_score)
    return target_probe_binds, non_target_probe_binds, target_primer_score, non_target_primer_score

shorter_short_list = []
targets_list: list = [file.stem.replace('.fna', '') for file in Path('targets').glob('*.fna.gz')]
non_targets_list: list = [file.stem.replace('.fna', '') for file in Path('non-targets').glob('*.fna.gz')]

for probescrape_line_dict in probescrape_shorter_list:

    target_probe_binds, non_target_probe_binds, target_primer_score, non_target_primer_score = _check_specificity(
        _probe_seq=probescrape_line_dict['probe_seq'],
        _forward_primer=probescrape_line_dict['fw_seq'],
        _reverse_primer=probescrape_line_dict['rev_seq'],
        _targets_list=targets_list,
        _non_targets_list=non_targets_list)

    probescrape_line_dict['probe_binds_ratio (target:non-target)'] = (target_probe_binds+1)/(non_target_probe_binds+1)
    
    probescrape_line_dict['target_probe_binds'] = target_probe_binds
    probescrape_line_dict['non_target_probe_binds'] = non_target_probe_binds

    probescrape_line_dict['target_primer_score'] = target_primer_score
    probescrape_line_dict['non_target_primer_score'] = non_target_primer_score

    shorter_short_list.append(probescrape_line_dict)
_output_results(shorter_short_list, 'probescrape-results.parsed.tested.csv')