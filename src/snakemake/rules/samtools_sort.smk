rule samtools_sort_extra:
    """
    Created:
        2017-05-06 22:50:39
    Aim:
        Sort bam, default by coordinates.
    Note:
        Sort bam by read name should be used for example before bedtools bamtobed bedpe.
    Test:
        out/samtools/sort/samtools/view_sam_to_bam/bowtie2/se_mm10/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR1202037.bam
        out/samtools/sort_-n/samtools/view_sam_to_bam/bowtie2/se_mm10/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR1202037.bam
    """
    input:
        bam="out/{filler}.bam"
    output:
        bam="out/{tool}{extra}/{filler}.bam"
    log:
        "out/{tool}{extra}/{filler}.log"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="samtools/sort"
    threads:
        MAX_THREADS
    conda:
        "../envs/samtools.yaml"
    envmodules:
        "samtools/1.14"
    shell:
        "samtools sort -@ {threads} -T {output.bam} {params.extra} {input.bam} -o {output.bam} &> {log}"

