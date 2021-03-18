rule gapfiller:
    """
    Created:
        2019-12-19 16:55:13
    Aim:
    Doc:
        https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.4.0/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php
    Note:
    Test:
        out/gatk4/HaplotypeCaller_fa-genome-GRCh38/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_se_fa-genome-GRCh38/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/Jurkat_SRR1057274_H3K27ac.vcf.gz
    """
    input:
        bam="out/{filler}.bam",
        fa= lambda wildcards: eval(mwconf['ids'][wildcards.fa_genome_id])
    output:
        vcf="out/{tool}{extra}_{fa_genome_id}/{filler}.vcf.gz"
    log:
        "out/{tool}{extra}_{fa_genome_id}/{filler}/log"
    benchmark:
        "out/{tool}{extra}_{fa_genome_id}/{filler}/benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "gapfiller/"
    threads:
        MAX_THREADS
    conda:
        "../envs/gapfiller.yaml"
    shell:
        "echo TODO"


