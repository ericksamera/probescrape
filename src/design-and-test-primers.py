from pathlib import Path
import subprocess
import csv
import random
from tracemalloc import start
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


def _parse_annotation(_target_region_str: str):
    """
    """
    chromosome_str, start_pos, end_pos = _parse_coords(_target_region_str)

    with open('probe_region', 'w') as probe_region_file: probe_region_file.write(f'{chromosome_str}\t{start_pos}\t{end_pos}')

    result = subprocess.run([
        'bedtools', 'intersect', '-a', 'reference/genomic.gff', '-b', 'probe_region'
    ], capture_output=True)
    result_lines = result.stdout.decode().split('\n')
    genes: list = []
    gene_products: list = []
    for line in list(filter(lambda line: 'product' in line, result_lines))[1:]:
        id_dict: dict = {item.strip().split('=')[0]: item.strip().split('=')[1] for item in line.strip().split('\t')[-1].split(';')}
        genes.append(id_dict.get('gene'))
        gene_products.append(id_dict.get('product'))

    print(set(genes), set(gene_products))
    return None

def main() -> None:

    with open('probescrape-results.csv') as input_csv:
        probe_results = []
        for line_dict in csv.DictReader(input_csv): 
            nt_flank: int = 50
            chromosome_str, start_pos, end_pos = _parse_coords(line_dict['target'])


            print(chromosome_str)
            _parse_annotation(chromosome_str, start_pos, end_pos)
            #target_region_str = f"{chromosome_str}:{start_pos-nt_flank}-{end_pos+nt_flank}"
            #_output_region_fasta(target_region_str)
            #_run_probesearch()
            #probe_results += _parse_probesearch_results(target_region_str, (end_pos+nt_flank)-(start_pos-nt_flank))
            #_output_results(probe_results); continue
        
        #_output_results(probe_results)

if __name__=='__main__':
    main()