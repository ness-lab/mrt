#!/usr/bin/env python3.6

import os


def create_output_directory(outpath):
    """Creates output directories, if they don't exist

    Args:
        outpath: Directory to which output will be written

    Returns:
        None: Creates directory, if necessary
    """

    if os.path.isdir(outpath):
        print("OUTPATH EXISTS, RUNNING SLiM SIMULATIONS")
    else:
        print("CREATING OUTPUT DIRECTORY")
        os.mkdir(outpath)
        print("RUNNING SLiM SIMULATIONS")
