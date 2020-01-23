rule bcftools_view_extin_to_extout:
    """
    Test:
        #out/bcftools/view/samtools/mpileup/samtools/index/samtools/sort/samtools/view_sam_to_bam/bowtie2/se_fasta_--rfg_1,1_-k_1_-q_-f_hg19/edena/assembling_-d_20_-c_20_-minCoverage_5/edena/overlapping/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_H3K27ac_contigs.vcf
        out/bcftools/view_vcf_to_bcf_-O_b/platypus/callVariants_--assemble=1_--minPosterior=0_fa-genome-hg19-main-chr_vcfgz-hg19-edena-contigs-with-inserts-gt-3bp-in-T11C-H3K27ac/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem2_se_bwa-index-hg19-main-chr-and-contigs-with-inserts-from-T11C-H3K27ac/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_H3K27ac.bcf
    """
    input:
        "out/{filler}.{extin}"
    output:
        "out/{tool}_{extin}_to_{extout}{extra}/{filler}.{extout}"
    params:
        extra = params_extra
    wildcard_constraints:
        extin="vcf|bcf",
        tool="bcftools/view"
    conda:
        "../../../../mw-lib/src/snakemake/envs/bcftools.yaml"
    shell:
        """
        bcftools view {params.extra} {input} > {output}
        """


