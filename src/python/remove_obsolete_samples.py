#!/usr/bin/env python

"""
Description:  This script look for files related to samples defined in a Sequencing_summary.xlsx and hardlink them into a tree structure.

Developer: Guillaume Charbonnier
Last modifications: 2023-02-12
Version: {v0.1}
"""

# import re  # Import regular expression
# import sys
import argparse
# import pysam
# import os.path
import os
import shutil
import glob
import pandas as pd
# from collections import defaultdict
__version__ = 0.1

def make_parser():
    """
    A function that returns the program parser.
    """

    parser = argparse.ArgumentParser(add_help=True)

    parser_grp_main = parser.add_argument_group('Arguments')

    parser_grp_main.add_argument

    parser_grp_main.add_argument(
        "-i",
        "--inp-dir",
        default = "out/ln/alias/sst/all_samples",
        help="The folder containing files to tidy."
    )

    parser_grp_main.add_argument(
        "-x",
        "--xlsx",
        type=str,
        help="The xlsx file containing the metadata to use to find samples and tidy them.",
        default="Sequencing_summary.xlsx",
        required=False)

    parser_grp_main.add_argument(
        "-b",
        "--by-column",
        nargs='+',
        type=str,
        help="The column names from the xlsx file to use to tidy.",
        default="sample_name",
        required=False)
        
    parser_grp_main.add_argument(
        "-d",
        "--delete",
        help="Delete file only this arg is used. Unsafe. Always run first without this argument and check all files listed to deletion.",
        default=False,
        type=bool,
    )

    return parser


def remove_obsolete_samples(
    inp_dir=None,
    xlsx=None,
    by_column=None,
    delete=None
):
    # Read the xlsx file with pandas
    samples = pd.read_excel(xlsx, sheet_name="samples", dtype=str)
    samples.fillna("NA", inplace=True)
    
    # Get the list of sample_name as file_prefixes
    file_prefixes = samples["sample_name"].tolist()

    # Get a list of all files in the current directory and subdirectories
    filepaths = []
    for root, dirs, files in os.walk(inp_dir):
        for file in files:
            filepaths.append(os.path.join(root, file))
    
    exception_suffixes = [
        "tsv",
        "html",
        "xlsx"
    ]
    # print(all_files)
    # Loop through the files and remove any that don't match a prefix
    for filepath in filepaths:
        if not matches_prefix(filepath, file_prefixes, exception_suffixes):
            print("file to remove : " + filepath)
            if delete:
                os.remove(filepath)
            pass



# Function to check if a file matches any of the prefixes
def matches_prefix(filepath, file_prefixes, exception_suffixes):
    filename = os.path.basename(filepath)
    for prefix in file_prefixes:
        # print("prefix : " + prefix)
        if filename.startswith(prefix):
            return True
    for suffix in exception_suffixes:
        if filename.endswith(suffix):
            return True
    return False


if __name__ == '__main__':

    args = make_parser()
    args = args.parse_args()
    args = dict(args.__dict__)
    remove_obsolete_samples(**args)
