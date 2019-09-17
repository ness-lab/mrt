#!/bin/sh
#SBATCH -N 1      # nodes requested
#SBATCH -n 1      # tasks requested
#SBATCH -c 9      # cores requested
#SBATCH --mem=400000  # memory in Mb
#SBATCH -o slim_sims_outfile  # send stdout to outfile
#SBATCH -e slim_sims_errfile  # send stderr to errfile

module load python/3.6.2
module load slim
module load gnu-parallel/20190222

parallel python ../python/SLiM_sims.py  ::: 100 1000 10000 ::: 0.1 0.5 1.0
