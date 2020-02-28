#!/usr/bin/env python3.6

import os


def create_output_directory(outpath):
    """Creates output directories, if they don't exist

    Args:
        outpath: Directory to which output will be written

    Returns:
        None: Creates directory, if necessary
    """

    print("Creating directory: {0}".format(outpath))
    if os.path.isdir(outpath):
        pass
    else:
        os.makedirs(outpath)
