#!/usr/bin/env python3.6

# Script uses subprocess to call SLiM on the command-line with varying
# values of N and the bottleneck proportion (as function of N). Allows
# simulation to easily be run on the cluster using GNU parallel.

import subprocess


def run_slim(N, bot):

    # Call SLiM from command line with N and bottleneck proportion values
    # (passed as command-line arguments)
    process = subprocess.Popen(["slim", "-d", "N=" + str(N), "-d", "bot=" + str(bot), "../slim/Bottleneck.slim"],
                               stdout = subprocess.PIPE,
                               stderr = subprocess.PIPE,
                               universal_newlines = True)
    out, err = process.communicate()

    # out is output of slim
    # err is error message from slim
    print(out)  # calls function to parse the output
    print(err)  # prints error message


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("N", help = "The desired population size", type = int)
    parser.add_argument("bot", help = "The desired strength of the population bottleneck. Expressed as the proportion of the population sampled during the bottleneck. 1.0 = No bottleneck", type = float)
    args = parser.parse_args()
    N = args.N
    bot = args.bot
    print(N, bot)

    run_slim(N, bot)
