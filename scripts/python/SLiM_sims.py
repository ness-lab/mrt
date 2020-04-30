#!/usr/bin/env python3.6

# Script uses subprocess to call SLiM on the command-line with varying
# values of N and the bottleneck proportion (as function of N). Allows
# simulation to easily be run on the cluster using GNU parallel.

import subprocess
import argparse
import os

from glob import glob
from tqdm import tqdm

from create_output_directory import create_output_directory
import slim_vcf2fasta_chrom


def args():
    """Parse and return command-line arguments

    Returns:
        Command-line arguments
    """
    parser = argparse.ArgumentParser(description='Run simulations in SLiM and convert VCFs to FASTA using Drosophila reference',
                                     usage='python3.7 SLiM_sims.py [options]')
    parser.add_argument('-n', '--pop_size', help='The desired population size', type=int, required=True)
    parser.add_argument('-b', '--bot', help='The desired strength of the population bottleneck. Expressed as the proportion of the population sampled during the bottleneck. 1.0=No bottleneck', type=float, required=True)
    parser.add_argument('-u', '--mutation_rate', help='Per base-pair mutation rate in scientific notation (e.g., 1e-8). (Default: 2.8e-9)', type=float, default=2.8e-9)
    parser.add_argument('-rec', '--recombination_rate', help='Per base-pair recombination rate rate in scientific notation (e.g., 1e-8). (Default: 2.14e-8)', type=float, default=2.14e-8)
    parser.add_argument('-s', '--slim_path', help='Path to SLiM script.', type=str, required=True)
    parser.add_argument('-f', '--fasta', required=True,
                        type=str, help='Reference sequence in FASTA format')
    parser.add_argument('-r', '--region', required=True,
                        type=str, help='samtools format region (1 index)')
    parser.add_argument('-m', '--mut_mat', required=True,
                        type=str, help='LDhelmet mut mat file')
    parser.add_argument('-o', '--outpath', help='Path to which VCFs from SLiM should be written', type=str, required=True)
    args = parser.parse_args()

    N = args.pop_size
    bot = args.bot
    outpath = str(args.outpath) + 'N{0}_bot{1}/'.format(N, bot)

    return N, bot, args.mutation_rate, args.recombination_rate, args.slim_path, args.fasta, args.region, args.mut_mat, outpath


def run_slim(N, bot, mu, rec, region, outpath, slim_path):
    """Run simulation using SLiM

    Args:
        N (int): Population size
        bot (float): Proportion of population to sample for bottleneck
        mu (float): Per-base-pair mutation rate
        rec (float): Per-base-pair recombination rate in represented as crossovers/bp (i.e., CM/bp x 10e-2)
        region (str): Genomic region in samtools format
        outpath (str): Path to which VCFs should be written
        slim_path (str): Path to SLiM script to execute at command-line

    Returns:
        None: VCFs for simulations written by SLiM script executed through command-line
    """
    # Calculate length of sequence from 'region' provided at command-line
    coords = region.split(':')[1]
    start, end = [int(n) for n in coords.split('-')]

    seq_length = (end - start) + 1

    print('Running SLiM simulations with N={0}, bot={1}, mu={2}, and rec={3}. Sequence length is {4}. VCFs in {5}'.format(N, bot, mu, rec, seq_length, outpath))

    seed = int(N) + int(bot * 100)

    # Call SLiM from command line with N and bottleneck proportion values
    # (passed as command-line arguments)
    outpath = "'" + outpath + "'"  # Required for command-line parsing and passing to SLiM
    process = subprocess.Popen(['slim', '-s', str(seed), '-d', 'N=' + str(N), '-d', 'bot=' + str(bot), '-d', 'mu=' + str(mu), '-d', 'rec=' + str(rec), '-d', 'seq_length=' + str(seq_length), '-d', 'outpath=' + str(outpath), slim_path],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               universal_newlines=True)
    out, err = process.communicate()

    # out is output of slim
    # err is error message from slim
    print(out)  # calls function to parse the output
    # print(err)  # prints error message


def find_vcfs(outpath, ext):
    """Run simulation using SLiM

    Args:
        ext (str): Extension of files to find (e.g., 'vcf')
        outpath (str): Path to where files are located

    Returns:
        process_find (obj of class: subprocess.Popen): Files on system with provided extension
    """
    print('Finding all {0} files in {1}'.format(ext, outpath))

    # Use find utility to identify all VCFs in outpath
    process_find = subprocess.Popen(['find', outpath, '-type', 'f',
                                    '-name', '*.' + ext],
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    universal_newlines=True)

    return process_find


def sort_vcfs(outpath):
    """Sort VCFs

    Args:
        outpath (str): Path to which sorted VCFs should be written

    Returns:
        None: Sorts VCFs in places through command-line
    """
    print('Sorting all VCFs in {0}'.format(outpath))

    total_vcfs = 0
    for vcf in tqdm(glob(outpath + '*.vcf')):

        filename = vcf.split('/')[-1].split('.vcf')[0]

        out_name = outpath + filename + '_sorted.vcf'

    # Files from 'find' are piped to 'xargs', which uses 'sort' to sort by position
        process_sort = subprocess.Popen(['vcf-sort', vcf],
                                        # stdin=process_find.stdout,
                                        stdout=open(out_name, 'w'),
                                        stderr=subprocess.PIPE,
                                        universal_newlines=True)

        out, err = process_sort.communicate()

        total_vcfs += 1

        # print(out)
        print(err)

        os.remove(vcf)

    print('Sorted a total of {0} VCFs'.format(total_vcfs))


def bgzip_vcfs(outpath):
    """BGzip VCFs

    Args:
        outpath (str): Path to which bgzipped VCFs should be written

    Returns:
        None: BGzipps VCFs in places through command-line
    """
    print('bgzipping all VCF files in {0}'.format(outpath))

    # Use find utility to identify all VCFs in outpath
    process_find = find_vcfs(outpath, 'vcf')

    # bgzip all found VCFs
    process_bgzip = subprocess.Popen(['xargs', '-n1', 'bgzip', '-f'],
                                     stdin=process_find.stdout,
                                     stderr=subprocess.PIPE,
                                     universal_newlines=True)

    out, err = process_bgzip.communicate()

    # print(out)
    print(err)


def tabix_vcfs(outpath):
    """Tabix index VCFs

    Args:
        outpath (str): Path to which tabixed VCFs should be written

    Returns:
        None: Tabix indexes VCFs in places through command-line
    """
    print('Tabix indexing all bgzipped VCF files in {0}'.format(outpath))

    # Use find utility to identify all VCFs in outpath
    process_find = find_vcfs(outpath, 'vcf.gz')

    # bgzip all found VCFs
    process_tabix = subprocess.Popen(['xargs', '-n1', 'tabix', '-f',
                                      '-p', 'vcf'],
                                     stdin=process_find.stdout,
                                     stderr=subprocess.PIPE,
                                     universal_newlines=True)

    out, err = process_tabix.communicate()

    # print(out)
    print(err)


def vcf2fasta(fasta, region, mut_mat, outpath):
    """Converts VCFs to FASTA

    Args:
        fasta (str): Path to reference genome in FASTA format
        region (str): Genomic region in samtools format
        mut_mat (str): Path to mutation matrix with substitution probabilities.
        outpath (str): Path to which FASTAs should be written

    Returns:
        None: Writes FASTA files to disk through command-line
    """
    fasta_outpath = outpath + 'fasta-files/'
    create_output_directory(fasta_outpath)

    total_vcfs = 0
    for vcf in tqdm(glob(outpath + '*.vcf.gz')):

        filename = fasta_outpath + vcf.split('/')[-1].split('.vcf.gz')[0] + '.fasta'

        slim_vcf2fasta_chrom.write_fasta(vcf, fasta, region, mut_mat, filename)

        total_vcfs += 1

    print('Converted a total of {0} VCF files from {1} to FASTA. FASTA files are stored in {2}'.format(total_vcfs, outpath, fasta_outpath))


if __name__ == '__main__':

    # Retrieve command-line arguments
    N, bot, mu, rec, slim_path, fasta, region, mut_mat, outpath = args()

    # Create output directory, if it doesn't exist
    create_output_directory(outpath)

    # Run simulations
    run_slim(N, bot, mu, rec, region, outpath, slim_path)

    # # Sort VCFs
    # sort_vcfs(outpath)

    # # bgzip files
    # bgzip_vcfs(outpath)

    # # Tabix index VCFs
    # tabix_vcfs(outpath)

    # # Convert VCFs to fasta
    # vcf2fasta(fasta, region, mut_mat, outpath)
