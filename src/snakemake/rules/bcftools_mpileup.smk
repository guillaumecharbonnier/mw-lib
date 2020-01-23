rule bcftools_mpileup_fagenome:
    """
    Note:
        Samtools mpileup is deprecated in favor of bcftools mpileup
    Aim:
        Intermediate step to produce VCF containing  TAL1 insertion:
        chr1·   47677744·   N·  2·  T+1Gt+21gggtaaaccgtctgttcagcg·  II

    Test:
        out/bcftools/mpileup_fa-genome-hg19-main-chr/samtools/index/samtools/sort/samtools/view_sam_to_bam/bowtie2/se_fasta_--rfg_1,1_-k_1_-q_-f_hg19/edena/assembling_-d_20_-c_20_-minCoverage_5/edena/overlapping/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_H3K27ac_contigs.vcf
        out/bcftools/mpileup_-I_fa-genome-hg19-main-chr/samtools/index/samtools/sort/samtools/view_sam_to_bam/bowtie2/se_fasta_--rfg_1,1_-k_1_-q_-f_hg19/edena/assembling_-d_20_-c_20_-minCoverage_5/edena/overlapping/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_H3K27ac_contigs.vcf

        out/bcftools/mpileup_-I_fa-genome-hg19-main-chr-and-contigs-with-insert-gt-3bp-in-T11C-H3K27ac/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem2_se_bwa-index-hg19-main-chr-and-contigs-with-inserts-gt-3bp-in-T11C-H3K27ac/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_H3K27ac.vcf
    """
    input:
        bam="out/{filler}.bam",
        fa= lambda wildcards: eval(config['ids'][wildcards.fa_genome_id])
    output:
        "out/{tool}{extra}_{fa_genome_id}/{filler}.vcf"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="bcftools/mpileup"
    conda:
        "../../../../mw-lib/src/snakemake/envs/bcftools.yaml"
    shell:
        """
        bcftools mpileup {params.extra} -f {input.fa} {input.bam} > {output}
        """

rule bcftools_mpileup_fagenome_bamlist:
    """
    Aim:
    Test:
        out/bcftools/mpileup_fa-genome-hg19-main-chr_bam-hg19-T11C-H3K27ac-dev1.vcf
    """
    input:
        bam = lambda wildcards: eval(config['ids'][wildcards.bam_list_id]),
        fa  = lambda wildcards: eval(config['ids'][wildcards.fa_genome_id])
    output:
        "out/{tool}{extra}_{fa_genome_id}_{bam_list_id}.vcf"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="bcftools/mpileup"
    conda:
        "../../../../mw-lib/src/snakemake/envs/bcftools.yaml"
    shell:
        """
        bcftools mpileup {params.extra} -f {input.fa} {input.bam} > {output}
        """






