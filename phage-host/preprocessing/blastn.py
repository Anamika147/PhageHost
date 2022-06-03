

import argparse
import multiprocessing
from pathlib import Path
import subprocess

desc = f'Run nucleotide BLAST for all files in a given directory'
p = argparse.ArgumentParser(description=desc)
p.add_argument('-i', '--in', metavar='indir', dest='indir', required=True
               )
p.add_argument('-d', '--db', metavar='dbdir', dest='dbdir', required=True
               )
p.add_argument('-o', '--out', metavar='outdir', dest='outdir', required=True
               )
p.add_argument('--task',  dest='task', choices=['blastn', 'megablast', 'dc-megablast'], 
               default='megablast')
p.add_argument('--num_threads', dest='num_threads', type=int,
               default=multiprocessing.cpu_count()
               )
args = p.parse_args()


in_path = Path(args.indir)
db_path = Path(args.dbdir).joinpath('db.fasta')
out_path = Path(args.outdir)
out_path.mkdir(parents=True, exist_ok=True)
for f in in_path.iterdir():
    genome_id = f.stem
    out_file = Path(out_path.joinpath(f'{genome_id}.txt'))
    cmd = [
        'blastn',
        '-task',
        args.task,
        '-query',
        f'{f}',
        '-db',
        f'{db_path}',
        '-outfmt',
        '6',
        '-out',
        f'{out_file}'
    ]
    subprocess.run(cmd)

