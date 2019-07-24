rule bbmap_reformat:
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
        tool = "bbmap/reformat"
    conda:
        "../envs/bbmap.yaml"
    shell:
        "reformat.sh in={input} out={output} {params.extra} &> {log}"

