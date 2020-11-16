rule samtools_sam_to_bam_bai_extra:
    """
    Created:
        2020-09-05 00:53:34
    Aim:
        Convert sam to bam. Additionnal filters can be applied using params_extra.
    Decoding SAM Flags:
        https://broadinstitute.github.io/picard/explain-flags.html
    Test:
        out/samtools/view_sam_to_bam_-q_30/bowtie2/se_mm10/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR1202037.bam
        out/samtools/view_sam_to_bam_-F_4_-q_5/bowtie2/se_mm10/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR1202037.bam
        out/samtools/view_sam_to_bam_-F_3588_-f_2_-q_20/
    """
    input:
        sam = "out/{filler}.sam"
    output:
        bam = "out/{tool}{extra}/{filler}.bam",
        bai = "out/{tool}{extra}/{filler}.bam.bai"
    log:
        "out/{tool}{extra}/{filler}.log"
    wildcard_constraints:
        tool="samtools/sam_to_bam_bai"
    params:
        extra = params_extra
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        (
        samtools view -bSh {params.extra} {input.sam} | samtools sort -o {output.bam}
        samtools index {output.bam}
        ) &> {log}
        """

