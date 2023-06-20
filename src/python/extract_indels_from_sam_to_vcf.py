#!/usr/bin/env python

import pysam
import argparse

parser = argparse.ArgumentParser()

parser.add_argument(
        "-s",
        "--sam",
        help="Input SAM")

parser.add_argument(
        "-g",
        "--genome",
        help="Input genome as fasta file")

parser.add_argument(
        "-v",
        "--vcf",
        help="Output VCF")


args = parser.parse_args()

# Load the reference genome sequence
reference = pysam.FastaFile(args.genome)

# Create a new VCF file
vcf_out = pysam.VariantFile(args.vcf, 'w')

# Add a header to the VCF file
vcf_out.header.add_line('##fileformat=VCFv4.2')
vcf_out.header.add_line('##INFO=<ID=CIGAR,Number=1,Type=String,Description="CIGAR string">')

samfile = pysam.AlignmentFile(args.sam, 'r')

# Add the chromosome names from the SAM file to the VCF header
for chrom in samfile.references:
    vcf_out.header.add_line(f'##contig=<ID={chrom}>')

for read in samfile.fetch():
    # Skip reads without a CIGAR string
    if read.cigartuples is None:
        continue

    # Skip reads without at least 20 matching bases at the beginning and end of the CIGAR string
    if read.cigartuples[0][0] != 0 or read.cigartuples[0][1] < 20 or read.cigartuples[-1][0] != 0 or read.cigartuples[-1][1] < 20:
        continue

    # Count the number of insertions and deletions in the CIGAR string
    insertion_count = 0
    deletion_count = 0
    for op, length in read.cigartuples:
        if op == 1:
            insertion_count += 1
        elif op == 2:
            deletion_count += 1

    # Skip reads with more than one insertion + deletion
    if insertion_count + deletion_count > 1:
        continue


    # Handle insertions and deletions by iterating through the CIGAR operations
    next_op_start = 0
    last_match_pos = None
    for op, length in read.cigartuples:
        if op == 0: # op == 0 corresponds to a match
            next_op_start += length
            last_match_pos = next_op_start - 1
        elif op == 1: # op == 1 corresponds to an insertion
            op_start = next_op_start
            op_end = next_op_start + length
            op_chromosome = read.reference_name
            op_seq = read.query_sequence[op_start:op_end]

            # Get the reference base immediately before the insertion
            ref_pos = read.reference_start + op_start - 1 # 0-based position
            ref_base = reference.fetch(op_chromosome, ref_pos, ref_pos + 1)

            # Create a new record for the VCF file
            record = vcf_out.new_record()
            record.chrom = op_chromosome
            record.pos = read.reference_start + op_start # VCF is 1-based
            record.id = '.'
            record.ref = ref_base
            record.alts = (ref_base + op_seq,)
            record.filter.add('PASS')

            # Add the CIGAR string to the record
            record.info['CIGAR'] = read.cigarstring

            # Write the record to the VCF file
            vcf_out.write(record)
        elif op == 2: # op == 2 corresponds to a deletion
            op_start = next_op_start
            op_end = next_op_start + length
            op_chromosome = read.reference_name

            # Get the reference sequence of the deletion and the base immediately before it
            ref_pos = read.reference_start + op_start - 1 # 0-based position
            ref_seq = reference.fetch(op_chromosome, ref_pos, ref_pos + length + 1)

            # Create a new record for the VCF file
            record = vcf_out.new_record()
            record.chrom = op_chromosome
            record.pos = read.reference_start + op_start # VCF is 1-based
            record.id = '.'
            record.ref = ref_seq
            record.alts = (ref_seq[0],)
            record.filter.add('PASS')

            # Add the CIGAR string to the record
            record.info['CIGAR'] = read.cigarstring

            # Write the record to the VCF file
            vcf_out.write(record)

# Close the VCF file
vcf_out.close()