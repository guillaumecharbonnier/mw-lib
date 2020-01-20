rule gatk4_haplotyper:
    """
    Created:
        2019-12-10 12:08:44
    Aim:
    Doc:
        https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.4.0/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php
    Note:
    Test:
        out/gatk4/HaplotypeCaller_fa-genome-GRCh38/samtools/index/picard/AddOrReplaceReadGroups_RGLB=rglbFiller_RGPL=illumina_RGPU=rgpuFiller_RGSM=rgsmFiller/samtools/sort/samtools/view_sam_to_bam/bwa/mem_se_fa-genome-GRCh38/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/Jurkat_SRR1057274_H3K27ac.vcf.gz

        out/gatk4/HaplotypeCaller_fa-genome-hg19-main-chr-and-contigs-with-insert-in-T11C-H3K27ac-v2/samtools/index/picard/AddOrReplaceReadGroups_RGLB\=rglbFiller_RGPL\=illumina_RGPU\=rgpuFiller_RGSM\=rgsmFiller/samtools/sort/samtools/view_sam_to_bam/bwa/mem2_se_bwa-index-hg19-main-chr-and-contigs-with-inserts-from-T11C-H3K27ac/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_H3K27ac.vcf.gz

        out/gatk4/HaplotypeCaller_-ERC_GVCF_fa-genome-hg19-main-chr-and-contigs-with-insert-in-T11C-H3K27ac-v2/samtools/index/picard/AddOrReplaceReadGroups_RGLB\=rglbFiller_RGPL\=illumina_RGPU\=rgpuFiller_RGSM\=rgsmFiller/samtools/sort/samtools/view_sam_to_bam/bwa/mem2_se_bwa-index-hg19-main-chr-and-contigs-with-inserts-from-T11C-H3K27ac/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_H3K27ac.vcf.gz

    """
    input:
        bam="out/{filler}.bam",
        fa= lambda wildcards: eval(config['ids'][wildcards.fa_genome_id])
    output:
        vcf="out/{tool}{extra}_{fa_genome_id}/{filler}.vcf.gz"
    log:
        "out/{tool}{extra}_{fa_genome_id}/{filler}/log"
    benchmark:
        "out/{tool}{extra}_{fa_genome_id}/{filler}/benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "gatk4/HaplotypeCaller"
    threads:
        16
    conda:
        "../envs/gatk4.yaml"
    shell:
        "gatk --java-options '-Xmx4g' HaplotypeCaller {params.extra} -I {input.bam} -R {input.fa} -O {output.vcf} &> {log}"


