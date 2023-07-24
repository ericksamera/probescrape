#!/usr/bin/env python3
__description__ =\
"""
Purpose: Given an input FASTA and region, add annotations from gff and vcf files to produce a GenBank file.
"""
__author__ = "Erick Samera"
__version__ = "1.0.0"
__comments__ = "stable enough"
# --------------------------------------------------
from argparse import (
    Namespace,
    ArgumentParser,
    RawTextHelpFormatter)
from pathlib import Path
# --------------------------------------------------
import primer3
import subprocess
import csv
import logging
import pandas as pd
from collections import Counter
# --------------------------------------------------
def get_args() -> Namespace:
    """ Get command-line arguments """

    parser = ArgumentParser(
        description=__description__,
        epilog=f"v{__version__} : {__author__} | {__comments__}",
        formatter_class=RawTextHelpFormatter)
    parser.add_argument('--targets',
        dest='targets_file_path',
        metavar="<INPUT TARGETS FILE>",
        type=Path,
        required=True,
        help="file containing target regions")
    parser.add_argument('--reference',
        dest='reference_fasta_path',
        metavar="<REFERENCE FASTA FILE>",
        type=Path,
        required=True,
        help="reference FASTA file")


    args = parser.parse_args()

    # parser errors and processing
    # --------------------------------------------------

    return args
# --------------------------------------------------
def _output_results(_results: list, _name_str: str) -> None: pd.DataFrame(_results).to_csv(f'{_name_str}', index=False); return None
def _output_region_fasta(_region_str: str, _reference_path: Path=Path('reference/GCA_016453205.2_ASM1645320v2_genomic.fna')) -> None:
    """
    """
    result = subprocess.run(['samtools', 'faidx', f'{_reference_path}', f'{_region_str}'], capture_output=True)
    parsed_seq: str = ''.join(result.stdout.decode().split('\n')[1:-1])
    return parsed_seq

def _get_products(_forward_primer_seq: str, _reverse_primer_seq: str, _mismatches: int, _targets_list: list, _non_targets_list: list) -> dict:
    """
    """
    with open('seq_primers', 'w') as primer_file: primer_file.write(f'000\t{_forward_primer_seq}\t{_reverse_primer_seq}\t50\t1000\n')

    result = subprocess.run([
        'ipcress',
        '-i', 'seq_primers',
        '--sequence', 'db_targets-and-non-targets.fna',
        '-m', f'{_mismatches}',
        '-p', 'FALSE',
        '-P', 'FALSE',
    ], capture_output=True)

    amplicons_list: list = ['_'.join(amplicon.split(':filter')[0].replace('ipcress: ', '').split('_')[:-1]) for amplicon in result.stdout.decode().split('\n')[:-2]]    
    #for i in amplicons_list: print(i)
    # print("All amplicons:", len(amplicons_list))
    # print("All amplicons, no duplicates:", len(Counter(amplicons_list).keys()))
    # print("Just M. bovis", len(list(filter(lambda x: 'bovis' in x, amplicons_list))))
    # print("Just M. bovis, no duplicates", len(Counter(filter(lambda x: 'bovis' in x, amplicons_list)).keys()))
    # print('---')
    # print(f'Mismatches: {_mismatches}')
    # print(f'{_forward_primer_seq}\t{_reverse_primer_seq}')

    ipcr_result: dict = {
        'all_amplicons_n': len(amplicons_list),
        'target_amplicons_n': len(list(filter(lambda x: 'bovis' in x, amplicons_list))),
        'amplicon_specific': len(amplicons_list)==len(list(filter(lambda x: 'bovis' in x, amplicons_list))),
        'amplicon_sensitive': len(Counter(filter(lambda x: 'bovis' in x, amplicons_list)).keys())/len(_targets_list)
    }

    return ipcr_result

def _design_primers(_input_seq: str) -> dict:
    """
    """
    primer_3_result = primer3.bindings.design_primers(
        seq_args={
            'SEQUENCE_ID': 'MH1000',
            'SEQUENCE_TEMPLATE': _input_seq,
            #'SEQUENCE_INCLUDED_REGION': [36,342]
        },
        global_args={
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_PICK_INTERNAL_OLIGO': 1,
            'PRIMER_INTERNAL_MAX_SELF_END': 8,
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 25,
            'PRIMER_MIN_TM': 57.0,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MAX_TM': 63.0,
            'PRIMER_MIN_GC': 20.0,
            'PRIMER_MAX_GC': 80.0,
            'PRIMER_MAX_POLY_X': 100,
            'PRIMER_INTERNAL_MAX_POLY_X': 100,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 50.0,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_MAX_SELF_ANY': 12,
            'PRIMER_MAX_SELF_END': 8,
            'PRIMER_PAIR_MAX_COMPL_ANY': 12,
            'PRIMER_PAIR_MAX_COMPL_END': 8,
            'PRIMER_PRODUCT_SIZE_RANGE': [[75, (len(_input_seq)-75)]],
        })

    return primer_3_result

def _parse_primer3_results(_primer_3_dict: dict) -> list:
    """
    """

    potential_primers: list = []
    results_num: list = [int(key.split('_')[2]) for key in _primer_3_dict if 'PRODUCT_SIZE' in key]
    for n in results_num:
        primer_result_dict = {
            'f_seq': _primer_3_dict[f'PRIMER_LEFT_{n}_SEQUENCE'],
            'f_tm': _primer_3_dict[f'PRIMER_LEFT_{n}_TM'],
            'r_seq': _primer_3_dict[f'PRIMER_RIGHT_{n}_SEQUENCE'],
            'r_tm': _primer_3_dict[f'PRIMER_RIGHT_{n}_TM'],
            'product_len': _primer_3_dict[f'PRIMER_PAIR_{n}_PRODUCT_SIZE'],
        }
        potential_primers.append(primer_result_dict)
    
    return potential_primers

def main() -> None:
    """
    """

    args = get_args()

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


    targets_list: list = [file.stem.replace('.fna', '') for file in Path('targets').glob('*.fna.gz')]
    non_targets_list: list = [file.stem.replace('.fna', '') for file in Path('non-targets').glob('*.fna.gz')]

    potential_primers: list = []

    potential_primer_results: list = []
    with open('seq-short-targets.csv', encoding='utf-8-sig') as seq_targets_csv:
        line_dicts = [line for line in csv.DictReader(seq_targets_csv)]
        for i_line_dict, line_dict in enumerate(line_dicts):
            parsed_seq: str = _output_region_fasta(line_dict['COORDS'], args.reference_fasta_path)
            if not len(parsed_seq) > 0: 
                logging.debug(f"Skipped region {line_dict['COORDS']} ({i_line_dict+1}/{len(line_dicts)}; {(i_line_dict+1)/(len(line_dicts)):.2f})).")
                continue
            primer_3_result: dict = _design_primers(parsed_seq)
            potential_primers += _parse_primer3_results(primer_3_result)
            logging.debug(f"Designed primers for region {line_dict['COORDS']} ({i_line_dict+1}/{len(line_dicts)}; {(i_line_dict+1)/(len(line_dicts)):.2f})).")
    
    logging.debug(f"Found {len(potential_primer_results)} potential primers.")
    for primer_pair in potential_primers:
        f_primer = primer_pair['f_seq']
        r_primer = primer_pair['r_seq']
        ipcr_results: dict = _get_products(
            f_primer,
            r_primer,
            _mismatches=2,
            _targets_list=targets_list,
            _non_targets_list=non_targets_list)
        primer_pair.update(ipcr_results)

        potential_primer_results.append(primer_pair)
        print(primer_pair)
    
    _output_results(potential_primer_results, 'seqscrape-results.csv')

#_parse_primer3_results(x)

"""
PRIMER_PAIR_0_PENALTY
PRIMER_LEFT_0_PENALTY
PRIMER_RIGHT_0_PENALTY
PRIMER_INTERNAL_0_PENALTY
PRIMER_LEFT_0_SEQUENCE
PRIMER_RIGHT_0_SEQUENCE
PRIMER_INTERNAL_0_SEQUENCE
PRIMER_LEFT_0
PRIMER_RIGHT_0
PRIMER_INTERNAL_0
PRIMER_LEFT_0_TM
PRIMER_RIGHT_0_TM
PRIMER_INTERNAL_0_TM
PRIMER_LEFT_0_GC_PERCENT
PRIMER_RIGHT_0_GC_PERCENT
PRIMER_INTERNAL_0_GC_PERCENT
PRIMER_LEFT_0_SELF_ANY_TH
PRIMER_RIGHT_0_SELF_ANY_TH
PRIMER_INTERNAL_0_SELF_ANY_TH
PRIMER_LEFT_0_SELF_END_TH
PRIMER_RIGHT_0_SELF_END_TH
PRIMER_INTERNAL_0_SELF_END_TH
PRIMER_LEFT_0_HAIRPIN_TH
PRIMER_RIGHT_0_HAIRPIN_TH
PRIMER_INTERNAL_0_HAIRPIN_TH
PRIMER_LEFT_0_END_STABILITY
PRIMER_RIGHT_0_END_STABILITY
PRIMER_PAIR_0_COMPL_ANY_TH
PRIMER_PAIR_0_COMPL_END_TH
PRIMER_PAIR_0_PRODUCT_SIZE
PRIMER_PAIR_0_PRODUCT_TM
"""

if __name__=="__main__":
    main()