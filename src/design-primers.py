from pathlib import Path
import subprocess

def _output_region_fasta(_region_str: str) -> None:
    """
    """
    result = subprocess.run([
        'samtools', 'faidx', 'reference/GCA_016453205.2_ASM1645320v2_genomic.fna', f'{_region_str}'
    ], capture_output=True)
    with open('target.fasta', mode='w') as fasta_output: fasta_output.write(result.stdout.decode())
    return None
def _run_primersearch(_probe_root: int, _probe_len: int) -> None:
    """
    """
    result = subprocess.run([
        'python', 'src-probedesign/primersearch.py', 
        'target.fasta',
        f'{_probe_root}',
        f'{_probe_len}',
        '--blastdb', 'db_targets-and-non-targets.fna',
        '--filter_tm',
        '--max_tm_diff-d', '5',
        '--min_primer_tm', '59',
        '--max_primer_tm', '64',
        '--mp_job', '12'
    ])
    return None

region: str ='CP058473.2:1069850-1070050'
probe_start: int = 78
probe_len: int = 18

_output_region_fasta(region)
_run_primersearch(probe_start, probe_len)