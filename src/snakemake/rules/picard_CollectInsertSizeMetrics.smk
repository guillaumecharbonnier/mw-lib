rule picard_CollectInsertSizeMetrics:
    """
    Created:
        2017-05-06 23:16:37
    Test:
        out/picard/CollectInsertSizeMetrics/samtools/index/samtools/sort/samtools/view_sam_to_bam/bowtie2/pe_hg19/sickle/pe_-t_sanger_-q_30/sra-tools/fastq-dump_pe/SRR1509753.pdf
    """
    input:
        bam="out/{filler}.bam"
    output:
        pdf="out/{tool}{extra}/{filler}.pdf",
        txt="out/{tool}{extra}/{filler}.txt"
    log:
            "out/{tool}{extra}/{filler}.log"
    benchmark:
            "out/{tool}{extra}/{filler}.benchmark.tsv"
    wildcard_constraints:
        tool="picard/CollectInsertSizeMetrics"
    params:
        extra=params_extra
    threads:
        1
    conda:
        "../envs/picard.yaml"
    shell:
        """
        picard CollectInsertSizeMetrics {params.extra} I={input.bam} O={output.txt} H={output.pdf} &> {log}
        """

