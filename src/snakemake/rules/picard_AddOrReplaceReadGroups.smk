rule picard_AddOrReplaceReadGroups:
    """
    Created:
        2017-05-06 23:16:37
    Doc:
        http://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups
    Test:
        out/picard/AddOrReplaceReadGroups_RGLB=rglbFiller_RGPL=illumina_RGPU=rgpuFiller_RGSM=rgsmFiller/samtools/index/samtools/sort/samtools/view_sam_to_bam/bowtie2/pe_hg19/sickle/pe_-t_sanger_-q_30/sra-tools/fastq-dump_pe/SRR1509753.bam
    """
    input:
        bam="out/{filler}.bam"
    output:
        bam="out/{tool}{extra}/{filler}.bam",
    log:
            "out/{tool}{extra}/{filler}.log"
    benchmark:
            "out/{tool}{extra}/{filler}.benchmark.tsv"
    wildcard_constraints:
        tool="picard/AddOrReplaceReadGroups"
    params:
        extra=params_extra
    threads:
        1
    conda:
        "../envs/picard.yaml"
    shell:
        """
        picard AddOrReplaceReadGroups {params.extra} I={input.bam} O={output.bam} &> {log}
        """
