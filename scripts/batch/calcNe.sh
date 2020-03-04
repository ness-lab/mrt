#!/bin/sh
#SBATCH -N 1      # nodes requested
#SBATCH -n 1      # tasks requested
#SBATCH -c 1      # cores requested
#SBATCH --mem=10000  # memory in Mb
#SBATCH -o calcNe_outfile  # send stdout to outfile
#SBATCH -e calcNe_errfile  # send stderr to errfile

module load python/3.6.2

python ../python/Ne_from_SLiM_vcfs.py -i ../../data/bottleneck/ -o ../../data/calcNe/
