#!/usr/bin/env python3
__description__ =\
"""
Purpose:
"""
__author__ = "Erick Samera"
__version__ = "0.0.1"
__comments__ = "stable enough"
# --------------------------------------------------
from argparse import (
    Namespace,
    ArgumentParser,
    RawTextHelpFormatter)
from pathlib import Path
# --------------------------------------------------
import subprocess
import Levenshtein
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
# --------------------------------------------------
def get_args() -> Namespace:
    """ Get command-line arguments """

    parser = ArgumentParser(
        description=__description__,
        epilog=f"v{__version__} : {__author__} | {__comments__}",
        formatter_class=RawTextHelpFormatter)
    parser.add_argument('-f', '--forward',
        dest='forward_seq',
        metavar="<FORWARD PRIMER>",
        type=str,
        required=True,
        help="forward primer sequence (5'-3')")
    parser.add_argument('-r', '--reverse',
        dest='reverse_seq',
        metavar="<REVERSE PRIMER>",
        type=str,
        required=True,
        help="reverse primer sequence (5'-3')")
    parser.add_argument('-m', '--mismatch',
        dest='mismatch_n',
        metavar="<MISMATCHES>",
        type=int,
        required=True,
        help="<n> mismatches to allow")

    args = parser.parse_args()

    # parser errors and processing
    # --------------------------------------------------

    return args
# --------------------------------------------------
def _get_products(_forward_seq: str, _reverse_seq: str, _mismatches: int, _targets_list: list, _non_targets_list: list) -> None:
    """
    """
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

    target_primer_score: float = 0
    non_target_primer_score: float = 0

    amplicons: list = []
    for seq_object in SeqIO.parse('products.fasta', 'fasta'):
        f_lev_distance: int = Levenshtein.distance(_forward_seq, seq_object.seq[:len(_reverse_seq)])
        r_lev_distance: int = Levenshtein.distance(Seq(_reverse_seq).reverse_complement(), seq_object.seq[-len(_forward_seq):])
        primer_score = ((len(_forward_seq)-f_lev_distance)+(len(_reverse_seq)-r_lev_distance))/(len(_forward_seq)+len(_reverse_seq))
        index_marker = 'GCA' if 'GCA' in seq_object.id else 'GCF'
        accession = '_'.join(seq_object.id.split('_')[:seq_object.id.split('_').index(index_marker)+2])
        amplicons.append(accession)
        if accession in _targets_list: target_primer_score += primer_score
        elif accession in _non_targets_list: non_target_primer_score += primer_score


    ipcr_results = {
        'represented_targets': len(set(list(filter(lambda accession: accession in _targets_list, amplicons)))),
        'represented_non_targets': len(set(list(filter(lambda accession: accession in _non_targets_list, amplicons)))),
        'total_amp_targets': len(list(filter(lambda accession: accession in _targets_list, amplicons))),
        'total_amp_non_targets': len(list(filter(lambda accession: accession in _non_targets_list, amplicons))),
        'target_primer_score': target_primer_score,
        'non_target_primer_score': non_target_primer_score,
        'target_primer_score_p': target_primer_score/len(_targets_list),
        'non_target_primer_score_p': non_target_primer_score/((len(_forward_seq)+len(_reverse_seq))*len(_non_targets_list))
    }

    for key, value in ipcr_results.items():
        print(f'{key}:\t{value}')

def _write_ipcr(_forward_seq: str, _reverse_seq: str, _min_product: int=50, _max_product: int=1000) -> None:
    """Function writes exonerate's ipcress primer file."""
    with open('primer', 'w') as primer_file: primer_file.write(f'M_BOVIS\t{_forward_seq}\t{_reverse_seq}\t{_min_product}\t{_max_product}\n')

def main() -> None:
    """
    """

    args = get_args()

    targets_list: list = [file.stem.replace('.fna', '') for file in Path('targets').glob('*.fna.gz')]
    non_targets_list: list = [file.stem.replace('.fna', '') for file in Path('non-targets').glob('*.fna.gz')]

    _write_ipcr(args.forward_seq, args.reverse_seq)
    _get_products(
        _forward_seq=args.forward_seq,
        _reverse_seq=args.reverse_seq,
        _mismatches=args.mismatch_n,
        _targets_list=targets_list,
        _non_targets_list=non_targets_list)

if __name__=='__main__':
    main()