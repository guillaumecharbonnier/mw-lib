rule snpeff:
    """
    Aim:
        First used to select subpopulations from bam.
        Trying this:
        https://www.biostars.org/p/225198/
    Note:
    Test:
        out/snpeff/GRCh38.99/platypus/callVariants_fa-genome-GRCh38/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_pe_fa-genome-GRCh38/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/803_H3K27ac.vcf
        out/snpeff/hg19/bcftools/setGT_somatic_to_missing/bcftools/fill-tags_-t_VAF/ln/updir/mw-tall-pediac-data/inp/vcf/3186patients_complet.sorted.vcf.gz
    """
    input:
        vcf="out/{filler}.{ext}",
    output:
        vcf="out/{tool}{extra}/{filler}.{ext}"
    log:
        "out/{tool}{extra}/{filler}.{ext}.log"
    benchmark:
        "out/{tool}{extra}/{filler}.{ext}.benchmark.tsv"
    params:
        extra = params_extra
    threads:
        MAX_THREADS
    wildcard_constraints:
        tool = "snpeff/",
        ext = "vcf|vcf.gz"
    conda:
        "../envs/snpeff.yaml"
    shell:
        """
        (
        if [ "{wildcards.ext}" == "vcf.gz" ]
        then
            zcat {input.vcf} | snpEff -Xmx8g {params.extra} | gzip -c > {output.vcf}
        fi
        if [ "{wildcards.ext}" == "vcf" ]
        then 
            snpEff -Xmx8g {params.extra} {input.vcf} > {output.vcf} 
        fi
        ) &> {log}
        """

rule snpsift_annotate_dbsnp:
    """
    Aim:
        Annotate a vcf with dbSnp. Notably, it will allow to identify novel from known variants in the vcf file.
    Test:
        out/snpsift/annotate_dbsnp/platypus/callVariants_fa-genome-GRCh38/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_pe_fa-genome-GRCh38/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/803_H3K27ac.vcf
        out/snpsift/annotate_dbsnp/python/extract_indels_from_sam_to_vcf_fa-genome-GRCh38-ensembl-r100/bowtie2/se_fa_Abraham2017_GRCh38-ensembl-r100/megahit/se/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac.contigs.vcf
    """
    input:
        vcf="out/{filler}.vcf",
        # Update using latest release on 2023-06-18
        # dbsnp="out/wget/https/ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz",
        # dbsnp_tbi="out/wget/https/ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.tbi"
        dbsnp = "out/bcftools/annnotate_rename_dbsnp_chrs_refseq_to_simple/wget/https/ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.vcf.gz",
        dbsnp_tbi = "out/bcftools/annnotate_rename_dbsnp_chrs_refseq_to_simple/wget/https/ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.vcf.gz.tbi"
        # dbsnp="out/wget/ftp/ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz",
        # dbsnp_tbi="out/wget/ftp/ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz.tbi"
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

rule snpsift_filter:
    """
    Aim:
        Initially created to filter vcf with no known id after annotating using dbsnp.
        Can be used for any filtering speicified in extra_ids.tsv
    Doc:
        http://pcingola.github.io/SnpEff/ss_filter/
    Test:
        out/snpsift/filter_no-id/snpsift/annotate_dbsnp/bcftools/merge/vcf-GRCh38-platypus-TALL-H3K27ac-H3K4me3.vcf
        out/snpsift/filter_high-moderate-impact/snpeff/hg19/bcftools/select_tall_samples/bcftools/setGT_somatic_to_missing/bcftools/fill-tags_-t_VAF/ln/updir/mw-tall-pediac-data/inp/vcf/3186patients_complet.sorted.vcf.gz
    """
    input:
        vcf="out/{filler}.{ext}"
    output:
        vcf="out/{tool}{extra}/{filler}.{ext}"
    log:
        "out/{tool}{extra}/{filler}.{ext}.log"
    benchmark:
        "out/{tool}{extra}/{filler}.{ext}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "snpsift/filter",
        ext = "vcf|vcf.gz"
    conda:
        "../envs/snpeff.yaml"
    shell:
        """
        (
        if [ "{wildcards.ext}" == "vcf.gz" ]
        then
            zcat {input.vcf} | SnpSift filter {params.extra} | gzip -c > {output.vcf}
        fi
        if [ "{wildcards.ext}" == "vcf" ]
        then 
            SnpSift filter {params.extra} {input.vcf} > {output.vcf}
        fi
        ) &> {log}
        """
