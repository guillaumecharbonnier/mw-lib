# import pysam

# def filter_reads(bamfile, snv_chrom, snv_pos, snv_allele):
#     filtered_reads = []
#     for read in bamfile.fetch(snv_chrom, snv_pos, snv_pos+1):
#         for pair in read.get_aligned_pairs():
#             read_pos, ref_pos = pair
#             if ref_pos == snv_pos and read.query_sequence[read_pos] == snv_allele:
#                 filtered_reads.append(read)
#     return filtered_reads

# bamfile = pysam.AlignmentFile("out/samtools/index/samtools/sort/samtools/view_sam_to_bam_-q_30/bowtie2/pe_GRCh38/sickle/pe_-t_sanger_-q_20/ln/alias/sst/all_samples/fastq/ATOJ7_CB548.bam", "rb")
# snv_chrom = "chr14"  # replace with your SNV chromosome
# snv_pos = 22448641  # replace with your SNV position (1-based)
# snv_pos = snv_pos - 1 # pysam uses 0-based coordinates
# snv_allele = "G"  # replace with your SNV allele

# filtered_reads = filter_reads(bamfile, snv_chrom, snv_pos, snv_allele)

# print(filtered_reads)



import pysam

def filter_reads(bamfile, snv_chrom, snv_pos, snv_allele):
    filtered_reads = set()
    # Adjust for 0-based coordinates used by pysam
    snv_pos_0based = snv_pos - 1
    for read in bamfile.fetch(snv_chrom, snv_pos_0based, snv_pos_0based + 1):
        for pair in read.get_aligned_pairs():
            read_pos, ref_pos = pair
            if ref_pos == snv_pos_0based and read.query_sequence[read_pos] == snv_allele:
                filtered_reads.add(read.query_name)
    return filtered_reads

bamfile = pysam.AlignmentFile("out/samtools/index/samtools/sort/samtools/view_sam_to_bam_-q_30/bowtie2/pe_GRCh38/sickle/pe_-t_sanger_-q_20/ln/alias/sst/all_samples/fastq/ATOJ7_CB548.bam", "rb")
snv_chrom = "chr14"  # replace with your SNV chromosome
snv_pos = 22448641  # replace with your SNV position (1-based)
snv_allele = "A"  # replace with your SNV allele

filtered_reads = filter_reads(bamfile, snv_chrom, snv_pos, snv_allele)

# Create output dir
import os
outdir = "out/py/filter_snv_allele"
os.makedirs(outdir, exist_ok=True)

outfilepath = os.path.join(outdir, "ATOJ7_CB548.bam")

# Write all reads except the filtered ones to a new BAM file
outfile = pysam.AlignmentFile(outfilepath, "wb", header=bamfile.header)
for read in bamfile.fetch():
    if read.query_name not in filtered_reads:
        outfile.write(read)
outfile.close()

# Index the new BAM file
pysam.index(outfilepath)