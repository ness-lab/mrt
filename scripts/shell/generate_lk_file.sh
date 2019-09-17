#!/bin/sh
#SBATCH -N 1      # nodes requested
#SBATCH -n 1      # tasks requested
#SBATCH -c 5     # cores requested
#SBATCH --mem=40000  # memory in Mb
#SBATCH -o generate_lk_outfile  # send stdout to outfile
#SBATCH -e generate_lk_errfile  # send stderr to errfile

~/Downloads/LDhat-master/complete -n 100 -rhomax 100 -n_pts 101 -theta $1 &&

cat new_lk.txt > $2/N100_lk_file.txt

rm new_lk.txt

