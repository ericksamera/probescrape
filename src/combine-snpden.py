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
import pandas as pd
# --------------------------------------------------
def get_args() -> Namespace:
    """ Get command-line arguments """

    parser = ArgumentParser(
        description=__description__,
        epilog=f"v{__version__} : {__author__} | {__comments__}",
        formatter_class=RawTextHelpFormatter)
    parser.add_argument('--targets',
        dest='targets_vcf',
        metavar="<TARGETS VCF>",
        type=Path,
        required=True,
        help="input FASTA file")
    parser.add_argument('--non-targets',
        dest='non_targets_vcf',
        metavar="<NON-TARGETS VCF>",
        type=Path,
        required=True,
        help="input FASTA file")
    parser.add_argument('--output',
        dest='output_file',
        metavar="<OUTPUT CSV FILE>",
        type=Path,
        required=True,
        help="output CSV file")
    parser.add_argument('--bin_size',
        dest='bin_size',
        metavar="<BIN SIZES>",
        type=int,
        required=True,
        default=100,
        help="bin size (bp)")

    args = parser.parse_args()

    # parser errors and processing
    # --------------------------------------------------
    if not Path(args.targets_vcf).is_file(): parser.error("asdasd")
    if not Path(args.non_targets_vcf).is_file(): parser.error("asdasd")

    return args
# --------------------------------------------------
def main() -> None:
    """ Main stuff. """
    
    args = get_args()

    snp_positions: dict = {}
    with open(args.targets_vcf) as targets_vcf_file:
        targets_vcf_file.readline()
        for line in targets_vcf_file.readlines():
            chromosome_str, bin_start, snp_count, _ = line.strip().split('\t')
            coords_str: str = f"{chromosome_str}:{int(bin_start)}-{int(bin_start)+args.bin_size}"
            if coords_str not in snp_positions: snp_positions[coords_str] = {
                'COORDS': coords_str,
                'TARGET': snp_count,
            }
    
    with open(args.non_targets_vcf) as targets_vcf_file:
        targets_vcf_file.readline()
        for line in targets_vcf_file.readlines():
            chromosome_str, bin_start, snp_count, _ = line.strip().split('\t')
            coords_str: str = f"{chromosome_str}:{int(bin_start)}-{int(bin_start)+args.bin_size}"
            if coords_str in snp_positions: snp_positions[coords_str]['NON-TARGET'] = snp_count

    pd.DataFrame.from_dict(snp_positions, orient='index').to_csv(args.output_file, index=False)
# --------------------------------------------------
if __name__ == '__main__':
    main()
