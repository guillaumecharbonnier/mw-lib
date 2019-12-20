rule platypus_callVariants:
    """
    Try ogap as alternative to bwa-mem

    Test:
        out/platypus/callVariants_fa-genome-GRCh38/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_se_fa-genome-GRCh38/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/Jurkat_SRR1057274_H3K27ac.vcf
        out/platypus/callVariants_fa-genome-hg19-main-chr/samtools/index/samtools/sort/samtools/view_sam_to_bam/bowtie2/se_fasta_--rfg_1,1_-k_1_-q_-f_hg19/edena/assembling_-d_20_-c_20_-minCoverage_5/edena/overlapping/gunzip/to-stdout/bowtie/se_--chunkmbs_256_--best_--strata_-m_1_-n_2_ebwt-hg19-main-chr/ln/alias/sst/all_samples/fastq/T11C_H3K27ac_unmapped_contigs.vcf

        bbmap/callvariants_fa-genome-GRCh38-r94-chr1/bbmap/se_fa-genome-GRCh38-r94-chr1/ln/alias/sst/all_samples/fastq/Jurkat_SRR1057274_H3K27ac.vcf

    """
    input:
        bam="out/{filler}.bam",
        bai="out/{filler}.bam.bai",
        fa= lambda wildcards: eval(config['ids'][wildcards.fa_genome_id])
    output:
        vcf="out/{tool}{extra}_{fa_genome_id}/{filler}.vcf"
    log:
        "out/{tool}{extra}_{fa_genome_id}/{filler}.log"
    benchmark:
        "out/{tool}{extra}_{fa_genome_id}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="platypus/callVariants"
    conda:
        "../envs/platypus.yaml"
    shell:
        """
        platypus callVariants --bamFiles={input.bam} --refFile={input.fa} --output={output.vcf} &> {log}
        """


