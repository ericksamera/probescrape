from pathlib import Path
import subprocess
import csv
import random
import pandas as pd
import logging
from Bio import SeqIO

def _output_region_fasta(_region_str: str) -> None:
    """
    """
    result = subprocess.run(['samtools', 'faidx', 'reference/GCA_016453205.2_ASM1645320v2_genomic.fna', f'{_region_str}'], capture_output=True)
    with open('target.fasta', mode='w') as fasta_output: fasta_output.write(result.stdout.decode())
    return None
def _parse_coords(_region_str: str) -> tuple:
    """
    """
    components = _region_str.split(':')
    chromosome_str = components[0]
    start_pos, end_pos = [int(num) for num in components[1].split('-')]
    return chromosome_str, start_pos, end_pos
def _run_probesearch(_region_str: str) -> None:
    """
    """
    logging.info(f'Running probesearch for region: {_region_str} .')
    result = subprocess.run([
        'python', 'src-probedesign/probesearch.py', 
        'target.fasta',
        #'--blastdb', 'db_targets-and-non-targets.fna',
        '--no_sens_spec_check',
        '--filter_tm',
        '--min_tm', '67',
        '--max_tm', '70',
        '--mp_job', '12'
    ])
    return None
def _mround(_number: int, _multiple: int) -> int: return _multiple * round(_number / _multiple)
def _flatten_list(_list: list) -> list: return [item for sublist in _list for item in sublist]
def _parse_probesearch_results(_target_region_str: str, _target_region_size: int, _region_buffer: int=50, _batch_k: int=3) -> None:
    """
    """
    left_region_buffer: int = _region_buffer
    right_region_buffer: int = _target_region_size - _region_buffer

    _batches_dict: dict = {}
    with open('probe_candidates.csv') as probe_csv:
        for line_dict in csv.DictReader(probe_csv):
            batch = _mround(int(line_dict['probe_root']), 25)
            if not (left_region_buffer < batch < right_region_buffer): continue
            if batch not in _batches_dict: _batches_dict[batch] = []
            _batches_dict[batch].append(line_dict)
    
    logging.info(f'Generated {len(_batches_dict)} batches of probes.')
    random_batch_choices = []
    logging.info(f'Picking {_batch_k} probes from each batch for a total of {len(_batches_dict)*_batch_k} probes.')
    for batch_num, batch in _batches_dict.items():
        random_batch_choices.append(random.choices(batch, k=_batch_k))
    probe_results = []
    for probe in _flatten_list(random_batch_choices):
        logging.info(f'Checking for potentially suitable primers.')
        _run_primersearch(probe['probe_root'], probe['probe_len'])
        if _check_primer_results(): 
            probe_to_add = probe
            probe_to_add['filesize'] = _check_primer_results()
            probe_to_add['target'] = _target_region_str
            probe_to_add['gene'], probe_to_add['gene_product'] = _find_nearby_annotations(_target_region_str)
            probe_results.append(probe_to_add)

    return probe_results
def _run_primersearch(_probe_root: int, _probe_len: int, _max_tm_diff: int=2, _min_primer_tm: int=59, _max_primer_tm: int=64, _mp_job: int=12) -> None:
    """
    """
    result = subprocess.run([
        'python', 'src-probedesign/primersearch.py', 
        'target.fasta',
        f'{_probe_root}',
        f'{_probe_len}',
        #'--blastdb', 'db_targets-and-non-targets.fna',
        '--no_sens_spec',
        '--filter_tm',
        '--max_tm_diff-d', f'{_max_tm_diff}',
        '--min_primer_tm', f'{_min_primer_tm}',
        '--max_primer_tm', f'{_max_primer_tm}',
        '--mp_job', f'{_mp_job}'
    ])
    return None
def _check_primer_results() -> None:
    """
    """
    if Path('primer_pairs.csv').stat().st_size > 145: return Path('primer_pairs.csv').stat().st_size
    else: return False
def _find_nearby_annotations(_target_region_str: str):
    """
    """
    chromosome_str, start_pos, end_pos = _parse_coords(_target_region_str)

    with open('probe_region', 'w') as probe_region_file: probe_region_file.write(f'{chromosome_str}\t{start_pos}\t{end_pos}')

    result = subprocess.run(['bedtools', 'intersect', '-a', 'reference/genomic.gff', '-b', 'probe_region'], capture_output=True)
    result_lines = result.stdout.decode().split('\n')
    genes: list = []
    gene_products: list = []
    for line in list(filter(lambda line: 'product' in line, result_lines))[1:]:
        id_dict: dict = {item.strip().split('=')[0]: item.strip().split('=')[1] for item in line.strip().split('\t')[-1].split(';')}
        genes.append(id_dict.get('gene', ''))
        gene_products.append(id_dict.get('product', ''))

    return ';'.join(sorted(set(genes))), ';'.join(sorted(set(gene_products)))
def _output_results(_results: list, _name_str: str) -> None: pd.DataFrame(_results).to_csv(f'{_name_str}', index=False); return None
def _generate_probes_short_list(_targets_file: str, _nt_flank: int=75, _batch_k: int=3) -> list:
    """
    """
    probe_results = []
    with open(_targets_file) as input_csv:
        for line_dict in csv.DictReader(input_csv): 
            chromosome_str, start_pos, end_pos = _parse_coords(line_dict['COORDS'])

            target_region_str = f"{chromosome_str}:{start_pos-_nt_flank if start_pos-_nt_flank > 0 else start_pos}-{end_pos+_nt_flank}"
            _output_region_fasta(target_region_str)
            _run_probesearch(target_region_str)
            probe_results += _parse_probesearch_results(target_region_str, (end_pos+_nt_flank)-(start_pos-_nt_flank), _batch_k=_batch_k)
        
    _output_results(probe_results, 'probescrape-results.csv')
    return probe_results
def _parse_probescrape_results(_probescrape_results: list, _batch_k: int=2, _primer_batch_mround: int=5) -> list:
    """
    """

    targets_list: list = [file.stem.replace('.fna', '') for file in Path('targets').glob('*.fna.gz')]
    non_targets_list: list = [file.stem.replace('.fna', '') for file in Path('non-targets').glob('*.fna.gz')]

    num_probes: int = len(_probescrape_results)
    logging.info(f'Found {num_probes} potential probes.')

    potential_primer_results: list = []
    for probescrape_line_dict in _probescrape_results:
        _output_region_fasta(probescrape_line_dict['target'])
        _run_primersearch(probescrape_line_dict['probe_root'], probescrape_line_dict['probe_len'])

        batches_dict: dict = {}
        with open('primer_pairs.csv') as probe_csv:
            for primer_line_dict in csv.DictReader(probe_csv):
                f_primer_batch = _mround(int(primer_line_dict['fw_root']), _primer_batch_mround)
                r_primer_batch = _mround(int(primer_line_dict['rev_root']), _primer_batch_mround)
                f_r_primer_batch = f"{f_primer_batch}-{r_primer_batch}"
                if f_r_primer_batch not in batches_dict: batches_dict[f_r_primer_batch] = []
                batches_dict[f_r_primer_batch].append(primer_line_dict)

        logging.info(f'Generated {len(batches_dict)} batches of primer pairs.')
        random_batch_choices = []
        logging.info(f'Picking {_batch_k} primer pairs from each batch for a total of {len(batches_dict)*_batch_k} primer pairs.')
        for batch_num, batch in batches_dict.items():
            random_batch_choices.append(random.choices(batch, k=_batch_k))
        
        for key in ('probe_root', 'probe_len', 'seq_rep', 'sens', 'spec', 'score', 'filesize'): probescrape_line_dict.pop(key)
        for primer_pair in _flatten_list(random_batch_choices):
            probescrape_line_dict['fw_seq'] = primer_pair['fw_seq']
            probescrape_line_dict['fw_tm'] = primer_pair['fw_tm']
            
            probescrape_line_dict['rev_seq'] = primer_pair['rev_seq']
            probescrape_line_dict['rev_tm'] = primer_pair['rev_tm']


            represented_targets, represented_non_targets, total_amp_targets, total_amp_non_targets = _get_products(
                _forward_primer=primer_pair['fw_seq'],
                _reverse_primer=primer_pair['rev_seq'],
                _mismatches=2,
                _targets_list=targets_list,
                _non_targets_list=non_targets_list)
            probescrape_line_dict['ratio (target:non-target)'] = (represented_targets+1)/(represented_non_targets+1)
            
            probescrape_line_dict['target_percentage_representation'] = represented_targets/len(targets_list)
            probescrape_line_dict['non_target_percentage_representation'] = represented_non_targets/len(non_targets_list)

            probescrape_line_dict['represented_targets'] = represented_targets
            probescrape_line_dict['represented_non_targets'] = represented_non_targets
            
            probescrape_line_dict['total_amplification_target'] = total_amp_targets
            probescrape_line_dict['total_amplification_non_target'] = total_amp_non_targets

            if probescrape_line_dict['target_percentage_representation'] < 0.99: continue
            potential_primer_results.append(probescrape_line_dict)
    
    _output_results(potential_primer_results, 'probescrape-results.parsed.csv')
    return potential_primer_results
def _get_products(_forward_primer: str, _reverse_primer: str, _mismatches: int, _targets_list: list, _non_targets_list: list) -> None:
    """
    """
    with open('primer', 'w') as primer_file: primer_file.write(f'M_BOVIS\t{_forward_primer}\t{_reverse_primer}\t50\t1000\n')

    result = subprocess.run([
        'ipcress',
        '-i', 'primer',
        '--sequence', 'db_targets-and-non-targets.fna',
        '-m', f'{_mismatches}',
        '-p', 'FALSE',
        '-P', 'FALSE',
    ], capture_output=True)

    amplicons_list: list = ['_'.join(amplicon.split(':filter')[0].replace('ipcress: ', '').split('_')[:-1]) for amplicon in result.stdout.decode().split('\n')[:-2]]
    represented_targets = len(set(list(filter(lambda accession: accession in _targets_list, amplicons_list))))
    represented_non_targets = len(set(list(filter(lambda accession: accession in _non_targets_list, amplicons_list))))
    total_amp_targets = len(list(filter(lambda accession: accession in _targets_list, amplicons_list)))
    total_amp_non_targets = len(list(filter(lambda accession: accession in _non_targets_list, amplicons_list)))
    return represented_targets, represented_non_targets, total_amp_targets, total_amp_non_targets
def _parse_broad_mismatch_probescrape_results(_probescrape_results: list) -> list:
    """
    """
    probescrape_results_non_redundant: dict = {}
    for line_dict in _probescrape_results:
        non_redundant_key = f"{line_dict['probe_seq']}_{line_dict['fw_seq']}_{line_dict['rev_seq']}"
        if non_redundant_key not in probescrape_results_non_redundant: probescrape_results_non_redundant[non_redundant_key] = line_dict
    probescrape_results_non_redundant: list = [result_dict for result_dict in probescrape_results_non_redundant.values()]

    probesearch_results_specific_primers: list = []
    targets_list: list = [file.stem.replace('.fna', '') for file in Path('targets').glob('*.fna.gz')]
    non_targets_list: list = [file.stem.replace('.fna', '') for file in Path('non-targets').glob('*.fna.gz')]
    for probescrape_line_dict in probescrape_results_non_redundant:
        for mismatch_num in (1, 0):
            represented_targets, represented_non_targets, total_amp_targets, total_amp_non_targets = _get_products(
                _forward_primer=probescrape_line_dict['fw_seq'],
                _reverse_primer=probescrape_line_dict['rev_seq'],
                _mismatches=3,
                _targets_list=targets_list,
                _non_targets_list=non_targets_list)
            print(represented_targets, represented_non_targets, total_amp_targets, total_amp_non_targets)
            probescrape_line_dict[f'{mismatch_num}_mismatch_ratio (target:non-target)'] = (represented_targets+1)/(represented_non_targets+1)
            
            probescrape_line_dict[f'{mismatch_num}_target_percentage_representation'] = represented_targets/len(targets_list)
            probescrape_line_dict[f'{mismatch_num}_non_target_percentage_representation'] = represented_non_targets/len(non_targets_list)

            probescrape_line_dict[f'{mismatch_num}_represented_targets'] = represented_targets
            probescrape_line_dict[f'{mismatch_num}_represented_non_targets'] = represented_non_targets
            
            probescrape_line_dict[f'{mismatch_num}_total_amplification_target'] = total_amp_targets
            probescrape_line_dict[f'{mismatch_num}_total_amplification_non_target'] = total_amp_non_targets

            if probescrape_line_dict[f'{mismatch_num}_target_percentage_representation'] < 0.99: continue
            probesearch_results_specific_primers.append(probescrape_line_dict)
    _output_results(probesearch_results_specific_primers, 'probescrape-results.parsed.primer-tested.csv')
    return probesearch_results_specific_primers
def _parse_primer_tested_results(_probescrape_results: list, _targets_list: list, _non_targets_list: list) -> list:
    """
    """
    _list_to_output: list = []
    for result_dict in _probescrape_results:

        with open('primer', 'w') as primer_file: primer_file.write(f"M_BOVIS\t{result_dict['fwd_seq']}\t{result_dict['rev_seq']}\t50\t1000\n")

        result = subprocess.run([
            'ipcress',
            '-i', 'primer',
            '--sequence', 'db_targets-and-non-targets.fna',
            '-m', '0',
            '-p', 'FALSE',
            '-P', 'TRUE',
            ], capture_output=True)

        with open('products.fasta', 'w') as products_file:
            for line in result.stdout.decode().split('\n'):
                if 'ipcress' in line: continue
                if line.startswith('>'): line = f">{line.strip().split('seq ')[1]}"
                products_file.write(line.strip()+'\n')
        
        target_binds: int = 0
        non_target_binds: int = 0
        for seq_object in SeqIO.parse('products.fasta', 'fasta'):
            if result_dict['probe_seq'] not in seq_object.seq: continue
            accession = '_'.join(seq_object.id.split('_')[:seq_object.id.split('_').index('GCA')+2])
            if accession in _targets_list: target_binds += 1
            elif accession in _non_targets_list: non_target_binds += 1
        
        result_dict['ratio_probe_binds'] = (target_binds+1)/(non_target_binds+1)
        result_dict['target_probe_binds'] = target_binds
        result_dict['non_target_probe_binds'] = non_target_binds
        _list_to_output.append(result_dict)
    _output_results(_list_to_output, 'probescrape-results.parsed.primer-tested.probe-tested.csv')
    return _list_to_output
def main() -> None:

    #args.output_path.joinpath('probesearch.log') = 
    logging.basicConfig(
        encoding='utf-8',
        level=logging.CRITICAL,
        handlers=[
            logging.FileHandler('probescrape.log'),
            logging.StreamHandler()
        ],
        format='%(asctime)s:%(levelname)s: %(message)s',
        datefmt='%m/%d/%Y %H:%M:%S',
    )

    probe_short_list: list = _generate_probes_short_list(_targets_file='', _batch_k=3)
    logging.critical(f'=== GENERATED A LIST OF POTENTIAL PROBE TARGETS ===')

    rough_tested_probes: list = _parse_probescrape_results(_probescrape_results=probe_short_list, _batch_k=3)
    logging.critical(f'=== GENERATED A LIST OF POTENTIAL PRIMERS FOR EACH PROBE TARGET ===')

    primer_tested_results: list = _parse_broad_mismatch_probescrape_results(_probescrape_results=rough_tested_probes)
    logging.critical(f'=== TESTED PRIMER SPECIFICITY FOR EACH PROBE TARGET ===')

    probe_tested_results: list = _parse_primer_tested_results(_probescrape_results=primer_tested_results)
    logging.critical(f'=== TESTED PROBE SPECIFICITY FOR EACH PROBE TARGET ===')
if __name__=='__main__':
    main()