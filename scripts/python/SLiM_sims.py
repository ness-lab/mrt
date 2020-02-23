#!/usr/bin/env python3.6

# Script uses subprocess to call SLiM on the command-line with varying
# values of N and the bottleneck proportion (as function of N). Allows
# simulation to easily be run on the cluster using GNU parallel.

import subprocess
import sys
import os


def create_output_directory(outpath):

    if os.path.isdir(outpath):
        print("OUTPATH EXISTS, RUNNING SLiM SIMULATIONS")
    else:
        print("CREATING OUTPUT DIRECTORY")
        os.mkdir(outpath)
        print("RUNNING SLiM SIMULATIONS")


def run_slim(N, bot, outpath, slim_path):

    # Call SLiM from command line with N and bottleneck proportion values
    # (passed as command-line arguments)
    outpath = "'" + outpath + "'"  # Required for command-line parsing and passing to SLiM
    process = subprocess.Popen(["slim", "-d", "N=" + str(N), "-d", "bot=" + str(bot), "-d", "outpath=" + str(outpath), slim_path],
                               stdout = subprocess.PIPE,
                               stderr = subprocess.PIPE,
                               universal_newlines = True)
    out, err = process.communicate()

    # out is output of slim
    # err is error message from slim
    # print(out)  # calls function to parse the output
    print(err)  # prints error message


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("N", help = "The desired population size", type = int)
    parser.add_argument("bot", help = "The desired strength of the population bottleneck. Expressed as the proportion of the population sampled during the bottleneck. 1.0 = No bottleneck", type = float)
    parser.add_argument("slim_path", help = "Path to SLiM script", type = str)
    parser.add_argument("outpath", help = "Path to which VCFs from SLiM should be written", type = str)
    args = parser.parse_args()

    # Retrieve command-line arguments
    N = args.N
    bot = args.bot
    slim_path = args.slim_path
    outpath = args.outpath

    # Change directory to script path so it can be run from any directory in project
    os.chdir(sys.path[0])

    # Create output directory, if it doesn't exit
    create_output_directory(outpath)

    # Run simulations
    run_slim(N, bot, outpath, slim_path)
