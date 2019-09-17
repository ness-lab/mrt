#!/bin/sh
#SBATCH -N 1      # nodes requested
#SBATCH -n 1      # tasks requested
#SBATCH -c 1     # cores requested
#SBATCH --mem=40000  # memory in Mb
#SBATCH -o create_csv_outfile  # send stdout to outfile
#SBATCH -e create_csv_errfile  # send stderr to errfile

module load python/2.7.14

python vcf_to_sfs_to_csv.py
