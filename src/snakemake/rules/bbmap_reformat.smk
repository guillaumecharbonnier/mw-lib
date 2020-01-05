rule bbmap_reformat_bam:
    """
    Created:
        2019-03-21 11:03:53
    Aim:
        First used to select subpopulations from bam.
    Note:
        Currently does not work because input bam is treated as unpaired...
    Test:
        out/bbmap/reformat_nuc-length/samtools/sort/samtools/merge_three_samples/samtools/merge_five_runs/samtools/view_sam_to_bam/bowtie2/pe_mm10/ln/pe_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run163_run167_run184_run187_run205/s1-MNase_Spm_WT_band1_s2-MNase_Spm_WT_band3_s3-MNase_Spm_WT_band6.bam
    """
    input:
        "out/{filler}"
    output:
        "out/{tool}{extra}/{filler}"
    log:
        "out/{tool}{extra}/{filler}.log"
    benchmark:
        "out/{tool}{extra}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "bbmap/reformat_bam"
    conda:
        "../envs/bbmap.yaml"
    shell:
        "reformat.sh in={input} out={output} {params.extra} &> {log}"

rule bbmap_reformat_se_extin_to_extout:
    """
    Created:
        2020-01-04 11:20:27
    Aim:
        Need to reformat this as fasta to use as reference input for GATK HaplotypeCaller:
        out/bowtie2/se_fasta_--rfg_1,1_-k_1_-q_-f_hg19/edena/assembling_-d_20_-c_20_-minCoverage_5/edena/overlapping/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_H3K27ac_contigs.sam
        But actually does not work as expected because it use only the first column as fasta sequence name whereas I would be more interested in naming using all columns excepted the sequence one. e.g.:
        out/edena/assembling_-d_20_-c_20_-minCoverage_5/edena/overlapping/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_H3K27ac_647339   16      chr20   33262100        255     42M     *       0       0       AGTGAGCCGAGATCGCGCCATTGCACTCCAGCCTGGGCAACA IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII      AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:42 YT:Z:UU
        TODO: Write awk script to do that.
        But first check if there are arguments to do that with bbmap reformat.
    Test:
        out/bbmap/reformat_se_sam_to_fasta/bowtie2/se_fasta_--rfg_1,1_-k_1_-q_-f_hg19/edena/assembling_-d_20_-c_20_-minCoverage_5/edena/overlapping/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_H3K27ac_contigs.fasta
    """
    input:
        "out/{filler}.{extin}"
    output:
        "out/{tool}_{extin}_to_{extout}{extra}/{filler}.{extout}"
    log:
        "out/{tool}_{extin}_to_{extout}{extra}/{filler}.log"
    benchmark:
        "out/{tool}_{extin}_to_{extout}{extra}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "bbmap/reformat_se",
        extin = "sam|fasta|fa", # Other ext can be added here if needed
        extout = "sam|fasta|fa",
    conda:
        "../envs/bbmap.yaml"
    shell:
        "reformat.sh in={input} out={output} {params.extra} &> {log}"

rule bbmap_reformat_fastq_pe:
    """
    Created:
        2019-12-20 15:35:11
    Aim:
        First use becaus abyss-pe want /1 and /2 suffixes on fastq
    Note:
    Test:
        out/bbmap/reformat_fastq_pe_addslash/ln/alias/sst/all_samples/fastq/489_H3K27ac_1.fastq.gz
    """
    input:
        m1="out/{filler}_1.{ext}",
        m2="out/{filler}_2.{ext}",
    output:
        m1="out/{tool}{extra}/{filler}_1.{ext}",
        m2="out/{tool}{extra}/{filler}_2.{ext}"
    log:
        "out/{tool}{extra}/{filler}.{ext}.log"
    benchmark:
        "out/{tool}{extra}/{filler}.{ext}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "bbmap/reformat_fastq_pe"
    conda:
        "../envs/bbmap.yaml"
    shell:
        "reformat.sh in={input.m1} in2={input.m2} out={output.m1} out2={output.m2} {params.extra} verifypaired=t &> {log}"

