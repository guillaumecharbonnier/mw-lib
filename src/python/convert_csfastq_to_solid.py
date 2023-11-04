#!/usr/bin/env python

from Bio import SeqIO
from Bio.SeqIO import FastaIO
import argparse

parser = argparse.ArgumentParser(description='Convert csfastq to native SOLID csfasta+qual.')
parser.add_argument('in_csfastq', type=str,
        help='Input csfastq file')
parser.add_argument('out_csfasta', type=str,
        help='Output csfasta file')
parser.add_argument('out_qual', type=str, help='Output qual file')

args = parser.parse_args()

fastq_sequences = SeqIO.parse(open(args.in_csfastq),'fastq-sanger')
with open(args.out_qual, "w") as out_qual:
    for fastq in fastq_sequences:
        qualities = str(fastq.format("qual"))
        # SOLID Native qual file have -1 instead of 0
        qualities = qualities.replace(' 0 ', ' -1 ')
        # Fastq file I have to process have this '!' not related to qual at the start of the qual string
        # '!' is converted to 0, hence this replacement
        qualities = qualities.replace('\n0 ', 'SUBSTITUTION_STRING')
        qualities = qualities.replace('\n', ' ')
        qualities = qualities.replace('SUBSTITUTION_STRING','\n')
        out_qual.write(qualities+"\n")

fastq_sequences = SeqIO.parse(open(args.in_csfastq),'fastq-sanger')
with open(args.out_csfasta, "w") as out_csfasta:
    fasta_out = FastaIO.FastaWriter(out_csfasta, wrap=None)
    fasta_out.write_file(fastq_sequences)
