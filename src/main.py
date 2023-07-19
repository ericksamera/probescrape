from pathlib import Path
import subprocess
import csv
import random
import pandas as pd

def _output_region_fasta(_region_str: str) -> None:
    """
    """
    result = subprocess.run([
        'samtools', 'faidx', 'reference/GCA_016453205.2_ASM1645320v2_genomic.fna', f'{_region_str}'
    ], capture_output=True)
    with open('target.fasta', mode='w') as fasta_output: fasta_output.write(result.stdout.decode())
    return None
def _parse_coords(_region_str: str) -> tuple:
    """
    """
    components = _region_str.split(':')
    chromosome_str = components[0]
    start_pos, end_pos = [int(num) for num in components[1].split('-')]
    return chromosome_str, start_pos, end_pos
def _run_probesearch() -> None:
    """
    """
    result = subprocess.run([
        'python', 'src-probedesign/probesearch.py', 
        'target.fasta',
        '--blastdb', 'db_targets-and-non-targets.fna',
        '--filter_tm',
        '--min_tm', '67',
        '--max_tm', '70',
        '--mp_job', '12'
    ])
    return None
def _mround(_number: int, _multiple: int) -> int: return _multiple * round(_number / _multiple)
def _flatten_list(_list: list) -> list: return [item for sublist in _list for item in sublist]
def _parse_probesearch_results(_target_region_str: str, _target_region_size: int, _region_buffer: int = 50) -> None:
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
    
    random_batch_choices = []
    for batch_num, batch in _batches_dict.items():
        random_batch_choices.append(random.choices(batch, k=2))
    
    probe_results = []
    for probe in _flatten_list(random_batch_choices):
        _run_primersearch(probe['probe_root'], probe['probe_len'])
        if _check_primer_results(): 
            probe_to_add = probe
            probe_to_add['filesize'] = _check_primer_results()
            probe_to_add['target'] = _target_region_str
            probe_results.append(probe_to_add)

    return probe_results
def _run_primersearch(_probe_root: int, _probe_len: int, _max_tm_diff: int=5, _min_primer_tm: int=59, _max_primer_tm: int=64, _mp_job: int=12) -> None:
    """
    """
    result = subprocess.run([
        'python', 'src-probedesign/primersearch.py', 
        'target.fasta',
        f'{_probe_root}',
        f'{_probe_len}',
        '--blastdb', 'db_targets-and-non-targets.fna',
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
def _output_results(_results: list) -> None:
    """
    """
    pd.DataFrame(_results).to_csv('probescrape-results.csv')
    return None

def main() -> None:
    with open('short-targets.csv') as input_csv:
        probe_results = []
        for line_dict in csv.DictReader(input_csv): 
            nt_flank: int = 50
            chromosome_str, start_pos, end_pos = _parse_coords(line_dict['COORDS'])

            try:
                target_region_str = f"{chromosome_str}:{start_pos-nt_flank}-{end_pos+nt_flank}"
                _output_region_fasta(target_region_str)
                _run_probesearch()
                probe_results += _parse_probesearch_results(target_region_str, (end_pos+nt_flank)-(start_pos-nt_flank))
            except: _output_results(probe_results); continue
        
        _output_results(probe_results)

if __name__=='__main__':
    main()