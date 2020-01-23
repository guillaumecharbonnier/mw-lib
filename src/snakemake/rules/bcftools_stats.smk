rule bcftools_stats_dev:
    input:
        "out/scalpel/discovery_--single_--bed_chr1:10000000-60000000_fa-genome-hg19-main-chr/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_se_fa-genome-hg19-main-chr/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/Jurkat_SRR1057274_H3K27ac/variants.indel.vcf",
        "out/scalpel/discovery_--single_--bed_chr1:10000000-60000000_fa-genome-hg19-main-chr/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_se_fa-genome-hg19-main-chr/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_H3K27ac/variants.indel.vcf",
        "out/scalpel/discovery_--single_--bed_chr1:10000000-60000000_fa-genome-hg19-main-chr/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_se_fa-genome-hg19-main-chr/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_all_DNA_samples/variants.indel.vcf"
    output:
        "out/bcftools/stats_dev/stats.out"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools stats {input} > {output}"

rule bcftools_plot_vcfstats:
    input:
        "out/bcftools/stats_dev/stats.out"
    output:
         touch("out/bcftools/plot-vcfstats/done")
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        plot-vcfstats -p out/bcftools/plot-vcfstats {input}
        """
