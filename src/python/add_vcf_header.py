#!/usr/bin/env python

import os
import argparse
import re 

# Define the header
header = """##fileformat=VCFv4.2
##source=FillerHeader
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO"""

# Create the parser
parser = argparse.ArgumentParser(description='Add a header to a VCF file if it is not present.')

# Add the arguments
parser.add_argument('input_vcf', type=str, help='The input VCF file')
parser.add_argument('output_vcf', type=str, help='The output VCF file')

# Parse the arguments
args = parser.parse_args()

# Regular expression for a simple VCF line check
vcf_line_regex = re.compile(r'^\S+\t\d+\t\S+\t[ACGTN]+\t[ACGTN,]+\t')

with open(args.input_vcf, 'r') as f_in, open(args.output_vcf, 'w') as f_out:
    lines = f_in.readlines()

    # Check if the file content looks like a VCF
    if not any(vcf_line_regex.match(line) for line in lines):
        raise ValueError(f"The content of the file {args.input_vcf} does not look like a VCF file.")

    # Check if the file already has a header
    if not lines[0].startswith('##fileformat'):
        # If not, add the header
        f_out.write(header + '\n')

    # Write the original lines to the output file
    f_out.writelines(lines)
