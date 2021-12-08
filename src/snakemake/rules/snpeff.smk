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
        "snpEff -Xmx8g {params.extra} {input.vcf} > {output.vcf} 2> {log}"

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
        tool = "snpsift/annotate_dbsnp",
    conda:
        "../envs/snpeff.yaml"
    shell:
        "SnpSift annotate {params.extra} {input.dbsnp} {input.vcf} > {output.vcf} 2> {log}"

rule wget_gwasCat:
    output:
        "out/wget_gwasCat/gwas_catalog_v1.0.2-associations_e104_r2021-11-22.tsv"
    shell:
        """
        cd out/wget_gwasCat
        wget 'https://www.ebi.ac.uk/gwas/api/search/downloads/alternative' --output-document 'gwas_catalog_v1.0.2-associations_e104_r2021-11-22.tsv'
        """

rule snpsift_gwasCat:
    """
    Aim:
        Annotate a vcf with dbSnp. Notably, it will allow to identify novel from known variants in the vcf file.
    Test:
        out/snpsift/gwasCat/platypus/callVariants_fa-genome-GRCh38/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_pe_fa-genome-GRCh38/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/803_H3K27ac.vcf
    """
    input:
        vcf="out/{filler}.vcf",
        gwasCat="out/wget_gwasCat/gwas_catalog_v1.0.2-associations_e104_r2021-11-22.tsv"
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
        tool = "snpsift/gwasCat",
    conda:
        "../envs/snpeff.yaml"
    shell:
        "SnpSift gwasCat {params.extra} {input.gwasCat} {input.vcf} > {output.vcf} 2> {log}"

rule snpsift_dbnsfp:
    """
    Aim:
        Annotate a vcf with dbnsfp.
    Test:
        out/snpsift/dbnsfp/platypus/callVariants_fa-genome-GRCh38/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_pe_fa-genome-GRCh38/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/803_H3K27ac.vcf
    """
    input:
        vcf="out/{filler}.vcf",
        dbnsfp="out/wget/https/snpeff.blob.core.windows.net/databases/dbs/GRCh38/dbNSFP_4.1a/dbNSFP4.1a.txt.gz",
        dbnsfp_tbi="out/wget/https/snpeff.blob.core.windows.net/databases/dbs/GRCh38/dbNSFP_4.1a/dbNSFP4.1a.txt.gz.tbi"
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
        tool = "snpsift/dbnsfp",
    conda:
        "../envs/snpeff.yaml"
    shell:
        "SnpSift dbnfsp {params.extra} -v -db {input.dbnsfp} {input.vcf} > {output.vcf} 2> {log}"
