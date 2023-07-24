from pathlib import Path
import subprocess
import csv
import random
import pandas as pd
import logging
import Levenshtein
from Bio import SeqIO
from Bio.Seq import Seq

# ===========================
def _mround(_number: int, _multiple: int) -> int: return _multiple * round(_number / _multiple)
def _flatten_list(_list: list) -> list: return [item for sublist in _list for item in sublist]
def _output_results(_results: list, _name_str: str) -> None: pd.DataFrame(_results).to_csv(f'{_name_str}', index=False); return None
# ===========================
def _parse_coords(_region_str: str) -> tuple:
    """
    """
    components = _region_str.split(':')
    chromosome_str = components[0]
    start_pos, end_pos = [int(num) for num in components[1].split('-')]
    return chromosome_str, start_pos, end_pos
def _output_region_fasta(_region_str: str, _reference_path: Path=Path('reference/GCA_016453205.2_ASM1645320v2_genomic.fna')) -> None:
    """
    """
    logging.debug(f'Wrote region {_region_str} to FASTA file.')
    result = subprocess.run(['samtools', 'faidx', f'{_reference_path}', f'{_region_str}'], capture_output=True)
    with open('target.fasta', mode='w') as fasta_output: fasta_output.write(result.stdout.decode())
    return None
# ===========================
def _run_probesearch(_region_str: str, _min_tm: int=67, _max_tm: int=70, _mp_job: int=12) -> None:
    """
    """
    result = subprocess.run([
        'python', 'src-probedesign/probesearch.py', 
        'target.fasta',
        #'--blastdb', 'db_targets-and-non-targets.fna',
        '--no_sens_spec_check',
        '--filter_tm',
        '--min_tm', f'{_min_tm}',
        '--max_tm', f'{_max_tm}',
        '--mp_job', f'{_mp_job}'
        ])
    return None
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
# ===========================
def _check_primer_results() -> None:
    """
    """
    if Path('primer_pairs.csv').stat().st_size > 145: return Path('primer_pairs.csv').stat().st_size
    else: return False
def _find_nearby_annotations(_target_region_str: str, _annotations_path: Path=Path('reference/genomic.gff')):
    """
    """
    chromosome_str, start_pos, end_pos = _parse_coords(_target_region_str)

    with open('probe_region', 'w') as probe_region_file: probe_region_file.write(f'{chromosome_str}\t{start_pos}\t{end_pos}')

    result = subprocess.run(['bedtools', 'intersect', '-a', f'{_annotations_path}', '-b', 'probe_region'], capture_output=True)
    result_lines = result.stdout.decode().split('\n')
    genes: list = []
    gene_products: list = []
    for line in list(filter(lambda line: 'product' in line, result_lines))[1:]:
        id_dict: dict = {item.strip().split('=')[0]: item.strip().split('=')[1] for item in line.strip().split('\t')[-1].split(';')}
        genes.append(id_dict.get('gene', ''))
        gene_products.append(id_dict.get('product', ''))

    return ';'.join(sorted(set(genes))), ';'.join(sorted(set(gene_products)))
def _parse_probesearch_results(_target_region_str: str, _target_region_size: int, _region_buffer: int=50, _batch_k: int=3) -> list:
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
        logging.debug(f'Checking for potentially suitable primers.')
        _run_primersearch(probe['probe_root'], probe['probe_len'])
        if _check_primer_results(): 
            probe_to_add = probe
            probe_to_add['target'] = _target_region_str
            probe_to_add['gene'], probe_to_add['gene_product'] = _find_nearby_annotations(_target_region_str)
            probe_results.append(probe_to_add)

    return probe_results

def _get_products(_probe_seq: str, _forward_primer_seq: str, _reverse_primer_seq: str, _mismatches: int, _targets_list: list, _non_targets_list: list) -> dict:
    """
    """
    with open('primer', 'w') as primer_file: primer_file.write(f'M_BOVIS\t{_forward_primer_seq}\t{_reverse_primer_seq}\t50\t1000\n')

    target_probe_binds: int = 0
    non_target_probe_binds: int = 0
    target_primer_score: float = 0
    non_target_primer_score: float = 0

    result = subprocess.run([
        'ipcress',
        '--input', 'primer',
        '--sequence', 'db_targets-and-non-targets.fna',
        '--mismatch', f'{_mismatches}',
        '--pretty', 'FALSE',
        '--products', 'TRUE',
    ], capture_output=True)
    with open('products.fasta', 'w') as products_file:
            for line in result.stdout.decode().split('\n'):
                if 'ipcress' in line: continue
                if line.startswith('>'): line = f">{line.strip().split('seq ')[1]}"
                products_file.write(line.strip()+'\n')

    amplicons_list: list = []
    for seq_object in SeqIO.parse('products.fasta', 'fasta'):
        f_lev_distance: int = Levenshtein.distance(_forward_primer_seq, seq_object.seq[:len(_forward_primer_seq)])
        r_lev_distance: int = Levenshtein.distance(Seq(_reverse_primer_seq).reverse_complement(), seq_object.seq[-len(_forward_primer_seq):])
        primer_score = (f_lev_distance+r_lev_distance)/(len(_forward_primer_seq)+len(_reverse_primer_seq))
        index_marker = 'GCA' if 'GCA' in seq_object.id else 'GCF'
        accession = '_'.join(seq_object.id.split('_')[:seq_object.id.split('_').index(index_marker)+2])
        amplicons_list.append(accession)
        if accession in _targets_list: 
            if _probe_seq in seq_object.seq: target_probe_binds += 1
            target_primer_score += primer_score
        elif accession in _non_targets_list:
            if _probe_seq in seq_object.seq: non_target_probe_binds += 1
            non_target_primer_score += primer_score

    represented_targets: int = len(set(list(filter(lambda accession: accession in _targets_list, amplicons_list))))
    represented_non_targets: int = len(set(list(filter(lambda accession: accession in _non_targets_list, amplicons_list))))
    total_amp_targets: int = len(list(filter(lambda accession: accession in _targets_list, amplicons_list)))
    total_amp_non_targets: int = len(list(filter(lambda accession: accession in _non_targets_list, amplicons_list)))

    products_dict: dict = {
        'ratio_target_non_target': (represented_targets+1)/(represented_non_targets+1),
        'represented_targets': represented_targets,
        'represented_targets_p': represented_targets/len(_targets_list),
        'represented_non_targets': represented_non_targets,
        'represented_non_targets_p': represented_non_targets/len(_non_targets_list),
        'total_amp_targets': total_amp_targets,
        'total_amp_non_targets': total_amp_non_targets,
        'probe_binds_target': target_probe_binds,
        'probe_binds_non_target': non_target_probe_binds,
        'primer_score_target': target_primer_score,
        'primer_score_non_target': non_target_primer_score
        }

    return products_dict

def _generate_probes_short_list(_targets_file: Path, _nt_flank: int=75, _batch_k: int=3) -> list:
    """
    """
    probe_results = []
    with open(_targets_file) as input_csv:
        for line_dict in csv.DictReader(input_csv): 
            chromosome_str, start_pos, end_pos = _parse_coords(line_dict['COORDS'])

            target_region_str = f"{chromosome_str}:{start_pos-_nt_flank if start_pos-_nt_flank > 0 else 1}-{end_pos+_nt_flank}"
            _output_region_fasta(target_region_str)
            
            logging.info(f'Running probesearch for region: {target_region_str} .')
            _run_probesearch(target_region_str)
            probe_results += _parse_probesearch_results(target_region_str, (end_pos+_nt_flank)-(start_pos-_nt_flank), _batch_k=_batch_k)
        
    _output_results(probe_results, 'probescrape-results.csv')
    return probe_results
def _parse_probescrape_results(_probescrape_results: list, _batch_k: int=2, _primer_batch_mround: int=5, _mismatches: int=1) -> list:
    """
    """

    targets_list: list = [file.stem.replace('.fna', '') for file in Path('targets').glob('*.fna.gz')]
    non_targets_list: list = [file.stem.replace('.fna', '') for file in Path('non-targets').glob('*.fna.gz')]

    num_probes: int = len(_probescrape_results)
    logging.info(f'Found {num_probes} potential probes.')

    filtered_probescrape_results: dict = {}
    for probescrape_line_dict in _probescrape_results:
        if probescrape_line_dict['probe_seq'] not in filtered_probescrape_results: filtered_probescrape_results[probescrape_line_dict['probe_seq']] = probescrape_line_dict
    filtered_probescrape_results: list = [item for item in filtered_probescrape_results.values()]

    potential_primer_results: list = []
    for result_i, probescrape_line_dict in enumerate(filtered_probescrape_results):
        logging.debug(f'Checking probe result {result_i+1}/{len(filtered_probescrape_results)} ({(result_i+1)/len(filtered_probescrape_results)}).')
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
        
        # remove old probesearch headers
        for key in ('probe_root', 'probe_len', 'seq_rep', 'sens', 'spec', 'score',): probescrape_line_dict.pop(key)

        primer_pairs_to_test: list = _flatten_list(random_batch_choices)
        for primer_pair_i, primer_pair in enumerate(primer_pairs_to_test):
            logging.debug(f'Checking probe result {primer_pair_i+1}/{len(primer_pairs_to_test)} ({(primer_pair_i+1)/len(primer_pairs_to_test)}).')
            probescrape_line_dict['fw_seq'] = primer_pair['fw_seq']
            probescrape_line_dict['fw_tm'] = primer_pair['fw_tm']
            
            probescrape_line_dict['rev_seq'] = primer_pair['rev_seq']
            probescrape_line_dict['rev_tm'] = primer_pair['rev_tm']

            logging.debug('Performing in-silico PCR.')
            products_dict = _get_products(
                _probe_seq=probescrape_line_dict['probe_seq'],
                _forward_primer_seq=primer_pair['fw_seq'],
                _reverse_primer_seq=primer_pair['rev_seq'],
                _mismatches=_mismatches,
                _targets_list=targets_list,
                _non_targets_list=non_targets_list)
            probescrape_line_dict.update(products_dict)

            if probescrape_line_dict['represented_targets_p'] < 0.99: continue

            potential_primer_results.append(probescrape_line_dict)
    
    _output_results(potential_primer_results, 'probescrape-results.parsed.csv')
    return potential_primer_results
# ===========================
def main() -> None:

    #args.output_path.joinpath('probesearch.log') = 
    logging.basicConfig(
        encoding='utf-8',
        level=logging.DEBUG,
        handlers=[
            logging.FileHandler('probescrape.log'),
            logging.StreamHandler()
        ],
        format='%(asctime)s:%(levelname)s: %(message)s',
        datefmt='%m/%d/%Y %H:%M:%S',
    )

    probe_short_list: list = _generate_probes_short_list(_targets_file=Path('super-short-targets.csv'), _batch_k=4)
    logging.critical(f'=== GENERATED A LIST OF POTENTIAL PROBE TARGETS ===')

    rough_tested_probes: list = _parse_probescrape_results(_probescrape_results=probe_short_list, _batch_k=4)
    logging.critical(f'=== GENERATED A LIST OF POTENTIAL PRIMERS FOR EACH PROBE TARGET ===')
if __name__=='__main__':
    main()