rule vcftoolz_compare_dev:
    """
    Created:
        2020-01-22 10:36:30
    Note:
        Untested rule because vcftoolz conda environment is not available right now.
        Alternative found for my needs using R VennDiagram and Upsetr.
    """
    input:
        "out/scalpel/discovery_--single_--bed_chr1:10000000-60000000_fa-genome-hg19-main-chr/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_se_fa-genome-hg19-main-chr/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/Jurkat_SRR1057274_H3K27ac/variants.indel.vcf",
        "out/scalpel/discovery_--single_--bed_chr1:10000000-60000000_fa-genome-hg19-main-chr/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_se_fa-genome-hg19-main-chr/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_H3K27ac/variants.indel.vcf",
        "out/scalpel/discovery_--single_--bed_chr1:10000000-60000000_fa-genome-hg19-main-chr/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_se_fa-genome-hg19-main-chr/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_all_DNA_samples/variants.indel.vcf"
    output:
        "out/vcftoolz/compare_dev/stats.out"
    conda:
        "../envs/vcftoolz.yaml"
    shell:
        "vcftoolz compare {input} > {output}"


