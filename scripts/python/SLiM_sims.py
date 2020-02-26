#!/usr/bin/env python3.6

# Script uses subprocess to call SLiM on the command-line with varying
# values of N and the bottleneck proportion (as function of N). Allows
# simulation to easily be run on the cluster using GNU parallel.

import subprocess
import argparse

from create_output_directory import create_output_directory


def run_slim(N, bot, outpath, slim_path):

    # Call SLiM from command line with N and bottleneck proportion values
    # (passed as command-line arguments)
    outpath = "'" + outpath + "'"  # Required for command-line parsing and passing to SLiM
    process = subprocess.Popen(["slim", "-d", "N=" + str(N), "-d", "bot=" + str(bot), "-d", "outpath=" + outpath, slim_path],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               universal_newlines=True)
    out, err = process.communicate()

    # out is output of slim
    # err is error message from slim
    # print(out)  # calls function to parse the output
    print(err)  # prints error message


def bgzip(outpath):

    # Find all VCFs in outpath
    process_find = subprocess.Popen(['find', outpath, '-name', '*.vcf'],
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    universal_newlines=True)

    # bgzip all found VCFs
    process_bgzip = subprocess.Popen(['xargs', '-n1', 'bgzip', '-f'],
                                     stdin=process_find.stdout,
                                     stderr=subprocess.PIPE,
                                     universal_newlines=True)

    out, err = process_bgzip.communicate()

    # print(out)
    print(err)


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

    # Create output directory, if it doesn't exit
    create_output_directory(outpath)

    # Run simulations
    run_slim(N, bot, outpath, slim_path)

    # bgzip files
    bgzip(outpath)
