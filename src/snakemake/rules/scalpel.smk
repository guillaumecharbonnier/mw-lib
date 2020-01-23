rule scalpel_discovery:
    """
    Created:
        2019-12-10 12:08:44
    Aim:
    Doc:
        https://sourceforge.net/p/scalpel/wiki/Manual/
    Note:
    Test:
        out/scalpel/discovery_--single_--bed_chr1:47185416-47335629_fa-genome-GRCh38/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_se_fa-genome-GRCh38/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/Jurkat_SRR1057274_H3K27ac/bamfile.bam
        out/scalpel/discovery_--single_--bed_chr1:47000000-48000000_fa-genome-hg19-main-chr/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_se_fa-genome-hg19-main-chr/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/Jurkat_SRR1057274_H3K27ac/bamfile.bam


    """
    input:
        bam="out/{filler}.bam",
        fa= lambda wildcards: eval(config['ids'][wildcards.fa_genome_id])
    output:
        "out/{tool}{extra}_{fa_genome_id}/{filler}/bamfile.bam"
    log:
        "out/{tool}{extra}_{fa_genome_id}/{filler}/log"
    benchmark:
        "out/{tool}{extra}_{fa_genome_id}/{filler}/benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "scalpel/discovery"
    threads:
        16
    conda:
        "../envs/scalpel.yaml"
    shell:
        "scalpel-discovery {params.extra} --numprocs {threads} --bam {input.bam} --ref {input.fa} --dir `dirname {output}` &> {log}"


rule scalpel_discovery_bed:
    """
    Created:
        2019-12-10 12:08:44
    Aim:
    Doc:
        https://sourceforge.net/p/scalpel/wiki/Manual/
    Note:
    Test:
        out/scalpel/discovery_--single_fa-genome-hg19-main-chr_bed-hg19-main-chr/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_se_fa-genome-hg19-main-chr/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/Jurkat_SRR1057274_H3K27ac/bamfile.bam out/scalpel/discovery_--single_fa-genome-hg19-main-chr_bed-hg19-main-chr/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_se_fa-genome-hg19-main-chr/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_H3K27ac/bamfile.bam out/scalpel/discovery_--single_fa-genome-hg19-main-chr_bed-hg19-main-chr/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_se_fa-genome-hg19-main-chr/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_all_DNA_samples/bamfile.bam
    """
    input:
        bam="out/{filler}.bam",
        bed= lambda wildcards: eval(config['ids'][wildcards.bed_id]),
        fa= lambda wildcards: eval(config['ids'][wildcards.fa_genome_id])
    output:
        "out/{tool}{extra}_{fa_genome_id}_{bed_id}/{filler}/bamfile.bam"
    log:
        "out/{tool}{extra}_{fa_genome_id}_{bed_id}/{filler}/log"
    benchmark:
        "out/{tool}{extra}_{fa_genome_id}_{bed_id}/{filler}/benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "scalpel/discovery"
    threads:
        16
    conda:
        "../envs/scalpel.yaml"
    shell:
        "scalpel-discovery {params.extra} --numprocs {threads} --bam {input.bam} --bed {input.bed} --ref {input.fa} --dir `dirname {output}` &> {log}"


