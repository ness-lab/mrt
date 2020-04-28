#!/bin/sh
#SBATCH -N 1      # nodes requested
#SBATCH -n 1      # tasks requested
#SBATCH -c 12      # cores requested
#SBATCH --mem=400000  # memory in Mb
#SBATCH -o bottleneck_outfile  # send stdout to outfile
#SBATCH -e bottleneck_errfile  # send stderr to errfile

module load python/3.6.2
module load slim
module load gnu-parallel/20190222
module load samtools/1.9
module load vcftools/0.1.16

parallel python ../python/SLiM_sims.py -n {1} -b {2} -s {3} -f {4} -r {5} -m {6} -o {7} \
::: 100 500 1000 5000 \
::: 0.1 0.5 1.0 \
::: ../slim/Bottleneck.slim \
::: ../../resources/dmel-reference-r6.32.fasta \
::: 2L:1000001-2000000 \
::: ../../resources/dmel_RAL_mut_mat.txt \
::: ../../data/bottleneck/
