rule bedtools_bamtobed_extra:
    """
    Created:
        2019-01-31 23:53:03
    Aim:
        Convert bam raw reads into bed raw reads for Dinup.
    Note:
        bai is not needed, moreover, this rule should be run on bam sorted by read name, thus preventing to produce a working bai as bai can only be produced on bam sorted by coordinates.
    Test:
        out/bedtools/bamtobed/samtools/sort_-n/samtools/view_sam_to_bam/bowtie2/se_mm10/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR1202037.bed
        out/bedtools/bamtobed_-tag_NM/samtools/sort_-n/samtools/view_sam_to_bam/bowtie2/se_mm10/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR1202037.bed
    """
    input:
        bam="out/{filler}.bam"
    output:
        "out/{tool}{extra}/{filler}.bed"
    log:
        "out/{tool}{extra}/{filler}.log"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="bedtools/bamtobed"
    conda:
        "../envs/bedtools.yaml"
    threads:
        1
    shell:
        "bedtools bamtobed {params.extra} -i {input.bam} > {output} 2> {log}"

