#!/usr/bin/env python3.6

# Script uses subprocess to call SLiM on the command-line with varying
# values of N and the bottleneck proportion (as function of N). Allows
# simulation to easily be run on the cluster using GNU parallel.

import subprocess
import argparse
import os

from glob import glob

from create_output_directory import create_output_directory


def run_slim(N, bot, outpath, slim_path):

    # Call SLiM from command line with N and bottleneck proportion values
    # (passed as command-line arguments)
    outpath = "'" + outpath + "'"  # Required for command-line parsing and passing to SLiM
    process = subprocess.Popen(["slim", "-s", "42", "-d", "N=" + str(N), "-d", "bot=" + str(bot), "-d", "outpath=" + str(outpath), slim_path],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               universal_newlines=True)
    out, err = process.communicate()

    # out is output of slim
    # err is error message from slim
    # print(out)  # calls function to parse the output
    print(err)  # prints error message


def find_vcfs(outpath):

    # Use find utility to identify all VCFs in outpath
    process_find = subprocess.Popen(['find', outpath, '-type', 'f',
                                     '-name', '*.vcf'],
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    universal_newlines=True)

    return(process_find)


def sort_vcfs(outpath):

    # Use find utility to identify all VCFs in outpath
    process_find = find_vcfs(outpath)

    # Files from 'find' are piped to 'xargs', which uses 'sort' to sort by position
    process_sort = subprocess.Popen(['xargs', '-I%', '-n1',
                                     'sort', '-n', '-k2',
                                     '%', '-o', '%'],
                                    stdin=process_find.stdout,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    universal_newlines=True)

    out, err = process_sort.communicate()

    # print(out)
    print(err)


def bgzip_vcfs(outpath):

    # Use find utility to identify all VCFs in outpath
    process_find = find_vcfs(outpath)

    # bgzip all found VCFs
    process_bgzip = subprocess.Popen(['xargs', '-n1', 'bgzip', '-f'],
                                     stdin=process_find.stdout,
                                     stderr=subprocess.PIPE,
                                     universal_newlines=True)

    out, err = process_bgzip.communicate()

    # print(out)
    print(err)


def convert_to_ldhat(outpath):

    # Path where LDhat input files will be stored
    ldhat_file_path = outpath + "ldhad_input_files/"

    # Create directory if it doesn't exist
    create_output_directory(ldhat_file_path)

    for vcf in glob(outpath + "*.vcf.gz"):

        filename = vcf.split('/')[-1].split('.vcf.gz')[0]
        # print(filename)

        process_vcftools = subprocess.Popen(['vcftools', '--gzvcf', vcf,
                                             '--chr', '1', '--ldhat',
                                             '--out', ldhat_file_path + filename],
                                            stdout=subprocess.PIPE,
                                            stderr=subprocess.PIPE,
                                            universal_newlines=True)

        out, err = process_vcftools.communicate()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("N", help="The desired population size", type=int)
    parser.add_argument("bot", help="The desired strength of the population bottleneck. Expressed as the proportion of the population sampled during the bottleneck. 1.0=No bottleneck", type=float)
    parser.add_argument("slim_path", help="Path to SLiM script.", type=str)
    parser.add_argument("outpath", help="Path to which VCFs from SLiM should be written", type=str)
    args = parser.parse_args()

    # Retrieve command-line arguments
    N = args.N
    bot = args.bot
    slim_path = args.slim_path
    outpath = str(args.outpath) + "N{0}_bot{1}/".format(N, bot)

    # Create output directory for VCFs, if it doesn't exit
    create_output_directory(outpath)

    # Run simulations
    run_slim(N, bot, outpath, slim_path)

    # Sort VCFs
    sort_vcfs(outpath)

    # bgzip files
    bgzip_vcfs(outpath)

    # Create LDhat '.sites' and '.locs' files
    convert_to_ldhat(outpath)
