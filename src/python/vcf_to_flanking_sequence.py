#!/usr/bin/env python

import pysam
import argparse

parser = argparse.ArgumentParser()

parser.add_argument(
        "-v",
        "--vcf",
        help="Input VCF")

parser.add_argument(
        "-g",
        "--genome",
        help="Input genome as fasta file")

parser.add_argument(
        "-f",
        "--flank",
        help="Number of nucleotide to extend around REF and ALT",
        default=20,
        type=int)

args = parser.parse_args()

# open vcf file
#vcf = pysam.VariantFile("input.vcf")
vcf = pysam.VariantFile(args.vcf)
#"out/MindTheGap/fill/MindTheGap/find_pe_fa-genome-GRCh38/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/856_H3K27ac.insertions.vcf")
# open fasta file
#genome = pysam.FastaFile("genome.fa")
genome = pysam.FastaFile(args.genome)
#"out/cat/assembly_ensembl/GRCh38.fa")

# define by how many bases the variant should be flanked
flank = args.flank

# iterate over each variant
for record in vcf:
    # extract sequence
    #
    # The start position is calculated by subtract the number of bases
    # given by 'flank' from the variant position. The position in the vcf file
    # is 1-based. pysam's fetch() expected 0-base coordinate. That's why we
    # need to subtract on more base.
    #
    # The end position is calculated by adding the number of bases
    # given by 'flank' to the variant position. We also need to add the length
    # of the REF value and subtract again 1 due to the 0-based/1-based thing.
    #
    # Now we have the complete sequence like this:
    # [number of bases given by flank]+REF+[number of bases given by flank]
    # Min and Max here to prevent requesting out of range coordinates after flanking
    # TODO: This Min Max likely
    #seq = genome.fetch(record.chrom,
    #        max(record.pos-1-flank, 1),
    #        min(record.pos-1+len(record.ref)+flank,
    #            genome.get_reference_length(record.chrom)))

    # Handling the special case where the indel is detected at the begining of a chr
    # Found for 525_H3K27ac_H3K4me3 for the first time
    # Avoid this error
    # ValueError: invalid coordinates: start (1) > stop (0)
    if record.pos == 1:
        seq_upstream = ""
    else:
        seq_upstream = genome.fetch(
            record.chrom,
            max(record.pos-1-flank, 1),
            record.pos-1
        )

    seq_downstream = genome.fetch(record.chrom,
            record.pos-1+len(record.ref),
            min(record.pos-1+len(record.ref)+flank,
                genome.get_reference_length(record.chrom)))

    seq = seq_upstream + record.ref + seq_downstream

    # Apply a filter to remove indel if it falls in homopolymer regions
    seq_lowercase_count=0
    seq_uppercase_count=0
    for i in seq:
        if(i.islower()):
            seq_lowercase_count=seq_lowercase_count+1
        elif(i.isupper()):
            seq_uppercase_count=seq_uppercase_count+1

    ref_lowercase_count=0
    ref_uppercase_count=0
    for i in record.ref:
        if(i.islower()):
            ref_lowercase_count=ref_lowercase_count+1
        elif(i.isupper()):
            ref_uppercase_count=ref_uppercase_count+1

    if seq_uppercase_count > seq_lowercase_count and ref_uppercase_count > ref_lowercase_count :
        # print out tab seperated columns:
        # CRHOM, POS, REF, ALT, flanking sequencing with variant given in the format '[REF/ALT]'
        #print('>', record.chrom, ':', record.pos, '_REF_', record.ref, '\n',
        #    seq[:flank], record.ref, seq[flank+len(record.ref):], '\n',
        #    '>', record.chrom, ':', record.pos, '_ALT_', record.alts[0], '\n',
        #    seq[:flank], record.alts[0], seq[flank+len(record.ref):],
        #    sep="")

        # It is required to trim seqname else Fimo will crash (likely when string is > 100 char)
        # Trimming to 80 because it is already long enough
        ref_seqname = '>' + str(record.chrom) + ':' + str(record.pos) + '_REF_' + record.ref
        ref_seqname = (ref_seqname[:79] + '..') if len(ref_seqname) > 80 else ref_seqname

        alt_seqname = '>' + str(record.chrom) + ':' + str(record.pos) + '_ALT_' + record.alts[0]
        alt_seqname = (alt_seqname[:79] + '..') if len(alt_seqname) > 80 else alt_seqname
        print(
            ref_seqname, '\n',
            seq_upstream, record.ref, seq_downstream, '\n',
            alt_seqname, '\n',
            seq_upstream, record.alts[0], seq_downstream,
            sep=""
        )
