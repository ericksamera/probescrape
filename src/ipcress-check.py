from pathlib import Path
import subprocess
from collections import Counter

forward_primer: str = 'ACGTGAGCCGCCAAGAACT'
reverse_primer: str = 'CAAGCCGCAGTTCAAGGAGAAA'

def _get_products() -> None:
    """
    """
    result = subprocess.run([
        'ipcress',
        '-i', 'primer',
        '--sequence', 'db_targets-and-non-targets.fna',
        '-m', '3',
        '-p', 'FALSE',
        '-P', 'TRUE',
    ], capture_output=True)

    with open('products.fasta', 'w') as products_file:
        for line in result.stdout.decode().split('\n'):
            if 'ipcress' in line: continue
            if line.startswith('>'): line = f">{line.strip().split('seq ')[1]}"
            products_file.write(line.strip()+'\n')
    
    print('Wrote products')

def run_ipcr() -> None:
    """
    """
    result = subprocess.run([
        'ipcress',
        '-i', 'primer',
        '--sequence', 'db_targets-and-non-targets.fna',
        '-m', '3',
        '-p', 'FALSE',
        '-P', 'FALSE',
    ], capture_output=True)

    amplicons_list: list = ['_'.join(amplicon.split(':filter')[0].replace('ipcress: ', '').split('_')[:-1]) for amplicon in result.stdout.decode().split('\n')[:-2]]    
    for i in amplicons_list: print(i)
    print("All amplicons:", len(amplicons_list))
    print("All amplicons, no duplicates:", len(Counter(amplicons_list).keys()))
    print("Just M. bovis", len(list(filter(lambda x: 'bovis' in x, amplicons_list))))
    print("Just M. bovis, no duplicates", len(Counter(filter(lambda x: 'bovis' in x, amplicons_list)).keys()))

with open('primer', 'w') as primer_file:
    primer_file.write(f'M_BOVIS\t{forward_primer}\t{reverse_primer}\t0\t1000\n')


run_ipcr()
_get_products()