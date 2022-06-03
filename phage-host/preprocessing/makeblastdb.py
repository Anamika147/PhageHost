

import argparse
from pathlib import Path
import subprocess as sub

desc = f'Create reference BLAST database from directory of FASTA sequences'
p = argparse.ArgumentParser(description=desc)
p.add_argument('-i', '--in', metavar='indir', dest='indir', required=True,
               )
p.add_argument('-o', '--out', metavar='outdir', dest='outdir', required=True,
               )
args = p.parse_args()

in_path = Path(args.indir)
out_path = Path(args.outdir)

# Create one FASTA file containing all the genomes from input directory.
out_path.mkdir(parents=True, exist_ok=True)
db_path = out_path.joinpath('db.fasta')
oh = open(db_path, 'w')
for f in in_path.iterdir():
    genome_id = f.stem
    n = 1
    fh = open(f)
    for line in fh:
        if line.startswith('>'):
            oh.write(f'>{genome_id}|{n}\n')
            n += 1
        else:
            oh.write(line)
    fh.close()
oh.close()

# Run makeblastdb.
#shell=True
commd = [
    'makeblastdb',
    '-in',
    f'{db_path}',
    '-dbtype',
    'nucl'
]
sub.run(commd)

# Remove db.fasta file
db_path.unlink()