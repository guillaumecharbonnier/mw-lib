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
        "-o",
        "--out-dir",
        help="The output folder where files will be put. Default 'inp-dir' to put new folders at the same level of the inp-dir.",
        default="inp-dir",
        required=False
    )
    # Replacing this arg with csv one
    #parser_grp_main.add_argument(
    #    "-x",
    #    "--xlsx",
    #    type=str,
    #    help="The xlsx file containing the metadata to use to find samples and tidy them.",
    #    default="Sequencing_summary.xlsx",
    #    required=False)

    parser_grp_main.add_argument(
        "-c",
        "--csv",
        type=str,
        help="The csv file containing the metadata to use to find samples and tidy them.",
        default="Sequency_summary.csv",
        required=False)


    parser_grp_main.add_argument(
        "-b",
        "--by-metadata",
        nargs='+',
        type=str,
        help="The column names from the xlsx file to use to tidy.",
        default=["exp", "chip_target","type","Sample_Project","investigator","cell_type","donor_id"],
        # "run" should not be added here becaused it is currently already added by the workflow
        required=False)
    
    parser_grp_main.add_argument(
        "-m",
        "--by-merge-metadata",
        nargs="+",
        type=str,
        default = None, 
        help='The columns names from the xlsx file to merge to use to tidy. Call one time for each needed merge. Example: -m "exp" "chip_target" -m "exp" "type"',
        action='append'
    )
    # parser_grp_main.add_argument(
    #     "-r",
    #     "--remove-out-dir",
    #     action="store_const",
    #     const=True,
    #     default=False,
    #     help="If set, remove the output folder before running. Ignore if out-dir is set to inp-dir",
    # )
    
    return parser

def tidy_samples(
    inp_dir=None,
    out_dir=None,
    #xlsx=None,
    csv=None,
    by_metadata=None,
    by_merge_metadata=None,
    remove_out_dir=None
):
    # Check if inp_dir equals "inp-dir"
    if out_dir == "inp-dir":
        # Set out_dir to the dirname of inp_dir
        out_dir = os.path.dirname(inp_dir)
    
    print(by_merge_metadata)

    # Commented because it is risky to remove out_dir like this
    # Someone could define a wrong out_dir and remove a lot of data
    # # Remove out_dir content if it exists
    # else if os.path.exists(out_dir) and remove_out_dir:
    #     print("debug")
    #     shutil.rmtree(out_dir)

    # Create out_dir if it does not exist
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Read the xlsx file with pandas
    #samples = pd.read_excel(xlsx, sheet_name="samples", dtype=str)
    samples = pd.read_csv(csv)

    # df = df.dropna(subset = "sample_name")
    samples.fillna("NA", inplace=True)

    # Remove the '.0' suffix from the 'run' column
    samples['run'] = samples['run'].astype(str).str.replace('.0', '', regex=False)

    # Iterate over by_merge_metadata to create the merged columns if not None
    if by_merge_metadata is not None:
        for merge_metadata in by_merge_metadata:
            # Combine merge_metadata strings into one
            merged_metadata_column = "_".join(merge_metadata)
            # Combine the columns into one
            samples[merged_metadata_column] = samples[merge_metadata].astype(str).apply('_'.join, axis=1)
            # Append the merged column name to by_metadata
            by_metadata.append(merged_metadata_column)
    
    # Currently bugged function
    # check_inclusion(samples)

    # # Get the list of samples
    # samples = samples["sample_name"].unique()

    # # Warn if there are duplicated samples
    # if len(samples) != len(samples["sample_name"]):
    #     print("Warning: there are duplicated samples in the xlsx file.")

    # Iterate over samples
    for index, row in samples.iterrows():
        tidy_sample(row, inp_dir, out_dir, by_metadata)

def tidy_sample(
    row,
    inp_dir,
    out_dir,
    by_metadata
):
    # Get working directory
    cwd = os.getcwd()

    search_dir = os.path.join(cwd, inp_dir)

    # Move to search dir
    os.chdir(search_dir)
    files = glob.glob("**/" + row['sample_name'] + ".*", recursive=True)
    os.chdir(cwd)

    if files == []:
        print("Warning: sample {} has no file.".format(row['sample_name']))
        return()
    
    # Iterate over files
    for file_path_suffix in files:
        # Define the input file path
        inp_file_path = os.path.join(cwd, inp_dir, file_path_suffix)
        # Iterate over by_metadata
        for metadata in by_metadata:
            # Split metadata value by comma
            metadata_values = row[metadata].split(",")
            print(metadata_values)
            # Iterate over metadata values
            for metadata_value in metadata_values:
                # Create the output file_path
                out_file_path = os.path.join(cwd, out_dir, "by_" + metadata, metadata_value, file_path_suffix)
                # Get the output file folder
                out_folder = os.path.dirname(out_file_path)
                # Create the output folder if it does not exist
                if not os.path.exists(out_folder):
                    os.makedirs(out_folder)
                
                # Check if file_path_suffix equal 'bw/quantile_normalized/879_H3K4me3_run281_over_TH91_SP8_H3K4me3.bw'
                # if file_path_suffix == 'bw/quantile_normalized/879_H3K4me3_run281_over_TH91_SP8_H3K4me3.bw':
                    # print("file_path_suffix: {}".format(file_path_suffix))

                # print(inp_file_path + " -> " + out_file_path)

                # Check if out_file_path exists:
                if os.path.exists(out_file_path):
                    print("Warning: out_file_path {} already exists. If you see this even with '--remove-out-dir', this is likely caused by different values of sample_name sharing a common regex. e.g. '879_H3K4me3' and '879_H3K4me3_run2'. It should not be an issue.".format(out_file_path))
                    continue

                # Create the hard link
                os.link(inp_file_path, out_file_path)


    # Get sample_name as a string
    # print(row['sample_name'])

    # if row['specie'] in ['human', 'Human', 'Homo_sapiens']:
    #     assemblies = ["GRCh38", "hg19"]
    # else:
    #     assemblies = ["mm10", "mm9"]
    
    # for assembly in assemblies:
    #     search_dir = os.path.join(cwd, inp_dir, assembly)

    #     # Check if search_dir exists
    #     if not os.path.exists(search_dir):
    #         # print("Warning: folder {} does not exist.".format(search_dir))
    #         return

    #     # Move to search dir
    #     os.chdir(search_dir)
    #     print(search_dir)

    #     # print(search_dir)
    #     # files = glob.glob(".*/" + row['sample_name'] + "_*")
    #     pattern = "**/" + row['sample_name'] + "_*"
    #     print(pattern)
    #     # files = glob.glob("**/730*", recursive=True)
    #     # print(files) 

    #     # # Print files if not empty
    #     # if files != []:
        

    #     # file_path = os.path.join(inp_dir, assembly, "bw", file_name)

    #     # # Check if file_path exists
    #     # if not os.path.exists(file_path):
    #     #     print("Warning: file {} does not exist.".format(file_path))
    #     #     return
        
    #     # for metadata in by_metadata:

    #     #     # Get the metadata value
    #     #     # print(metadata)
    #     #     metadata_value = row[metadata]
    #     #     # Create the output folder
    #     #     out_folder = os.path.join(out_dir, assembly, metadata, metadata_value, "bw")
    #     #     print(out_folder)
    #     #     if not os.path.exists(out_folder):
    #     #         os.makedirs(out_folder)
    #     #     # Create the output file path
    #     #     out_file_path = os.path.join(out_folder, file_name)
    #     #     # Create the hard link
    #     #     os.link(file_path, out_file_path)
        

def check_inclusion(lst):
    for i in range(len(lst)):
        for j in range(i + 1, len(lst)):
            if isinstance(lst[i], str) and isinstance(lst[j], str) and lst[i] in lst[j]:
                return True
    return False


if __name__ == '__main__':

    args = make_parser()
    args = args.parse_args()
    args = dict(args.__dict__)
    tidy_samples(**args)
