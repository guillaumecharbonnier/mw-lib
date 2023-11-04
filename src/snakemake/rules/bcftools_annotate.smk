rule bcftools_annnotate_rename_dbsnp_chrs_refseq_to_simple:
    """
    Aim:
        I wanted to update the neoenhancer workflow to use the latest dbsnp vcf. However, they switched from using 1,2,3 chromosomes to NC000001, etc...
        I need to rename the chromosomes in the vcf.
        Note I wrote a very specific rule because I also need to append the vcf suffix to the output file.
    Test:
        out/bcftools/annotate_rename_dbsnp_chrs_refseq_to_simple/wget/https/ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.vcf.gz
    """
    input:
        vcf  = "out/wget/https/ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz",
        chr_map = "../mw-lib/src/tsv/GCF_000001405.40_primary_assembly_chromosomes.tsv"
    output:
        vcf = "out/bcftools/annnotate_rename_dbsnp_chrs_refseq_to_simple/wget/https/ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.vcf.gz",
        tbi = "out/bcftools/annnotate_rename_dbsnp_chrs_refseq_to_simple/wget/https/ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.vcf.gz.tbi"
    log:
        "out/bcftools/annnotate_rename_dbsnp_chrs_refseq_to_simple/wget/https/ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.log"
    benchmark:
        "out/bcftools/annnotate_rename_dbsnp_chrs_refseq_to_simple/wget/https/ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.benchmark.tsv"
    conda:
        "../../../../mw-lib/src/snakemake/envs/bcftools_tabix.yaml"
    shell:
        """
        bcftools annotate --rename-chrs {input.chr_map} {input.vcf} -Oz -o {output.vcf} &> {log}
        tabix -p vcf {output.vcf}
        """

