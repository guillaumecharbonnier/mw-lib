rule snpeff:
    """
    Aim:
        First used to select subpopulations from bam.
        Trying this:
        https://www.biostars.org/p/225198/
    Note:
    Test:
        out/snpeff/GRCh38.99/platypus/callVariants_fa-genome-GRCh38/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_pe_fa-genome-GRCh38/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/803_H3K27ac.vcf
    """
    input:
        vcf="out/{filler}.vcf",
    output:
        vcf="out/{tool}{extra}/{filler}.vcf"
    log:
        "out/{tool}{extra}/{filler}.log"
    benchmark:
        "out/{tool}{extra}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    threads:
        MAX_THREADS
    wildcard_constraints:
        tool = "snpeff/",
    conda:
        "../envs/snpeff.yaml"
    shell:
        "snpEff -Xmx8g {params.extra} {input.vcf} > {output.vcf} &> {log}"

rule snpsift_annotate_dbsnp:
    """
    Aim:
        Annotate a vcf with dbSnp. Notably, it will allow to identify novel from known variants in the vcf file.
    Test:
        out/snpsift/annotate/platypus/callVariants_fa-genome-GRCh38/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_pe_fa-genome-GRCh38/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/803_H3K27ac.vcf
    """
    input:
        vcf="out/{filler}.vcf",
        dbsnp="out/wget/ftp/ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz",
        dbsnp_tbi="out/wget/ftp/ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz.tbi"
    output:
        vcf="out/{tool}{extra}/{filler}.vcf"
    log:
        "out/{tool}{extra}/{filler}.log"
    benchmark:
        "out/{tool}{extra}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    threads:
        MAX_THREADS
    wildcard_constraints:
        tool = "snpsift/annotate",
    conda:
        "../envs/snpeff.yaml"
    shell:
        "SnpSift annotate {params.extra} {input.dbsnp} {input.vcf} > {output.vcf} &> {log}"
