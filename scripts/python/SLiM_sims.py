#!/usr/bin/env python3.6

# Script uses subprocess to call SLiM on the command-line with varying
# values of N and the bottleneck proportion (as function of N). Allows
# simulation to easily be run on the cluster using GNU parallel.

import subprocess
import sys


def parse(out):
    lines = out.split("\n")
    for line in lines:
        print(line)
        print("--------------")


def main():

    N = int(sys.argv[1])  # initial N value
    bottleneck = float(sys.argv[2])  # proportion of N

    # Call SLiM from command line with N and bottleneck proportion values
    # (passed as command-line arguments)
    process = subprocess.Popen(["slim", "-d", "N=" + str(N), "-d", "bottleneck=" + str(bottleneck), "../slim/Bottleneck.slim"],
                               stdout = subprocess.PIPE,
                               stderr = subprocess.PIPE,
                               universal_newlines = True)
    out, err = process.communicate()

    # out is output of slim
    # err is error message from slim
    parse(out)  # calls function to parse the output
    print(err)  # prints error message


if __name__ == "__main__":
    main()
