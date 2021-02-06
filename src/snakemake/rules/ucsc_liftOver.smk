rule ucsc_liftOver:
    """
    Created:
        2019-02-18 17:00:18
    Aim:
        Converting wig from danpos wiq to bigwig.
    Test:
        out/ucsc/liftOver_chain-mm10-to-mm9/cut/_-f1-6/macs2/callpeak_--broad/samtools/index/samtools/sort/samtools/view_sam_to_bam_-q_30/bowtie2/se_mm10/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR3126243_over_SRR3126242_peaks.bed
    """
    input:
        bed = "out/{filler}.bed",
        chain = lambda wildcards: eval(mwconf['ids'][wildcards.chain_id])
    output:
        bed = "out/ucsc/liftOver_{chain_id}/{filler}.bed",
        unmapped = "out/ucsc/liftOver_{chain_id}/{filler}.bed.unmapped"
    conda:
        "../envs/ucsc_liftover.yaml"
    shell:
        "liftOver {input.bed} {input.chain} {output.bed} {output.unmapped}"
