#!/bin/sh
#SBATCH -N 1      # nodes requested
#SBATCH -n 1      # tasks requested
#SBATCH -c 30      # cores requested
#SBATCH --mem=100000  # memory in Mb
#SBATCH -o ldhelmet_outfile  # send stdout to outfile
#SBATCH -e ldhelmet_errfile  # send stderr to errfile

module load ldhelmet/1.10
module load boost/1.6.5
module load gsl/2.6

sh ../shell/run-ldhelmet.sh ../../data/bottleneck/ ../../data/bottleneck/ldhelmet
