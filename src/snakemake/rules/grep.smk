rule grep_extra:
    """
    Created:
        2019-02-01 11:50:39
    Aim:
    Test:
        out/grep/ENSG00000277734/sed/add_chr/gunzip/to-stdout/wget/ftp/ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.gtf
        out/pbgzip_tabix/grep/INDEL/bcftools/mpileup_fa-genome-hg19-main-chr/samtools/index/samtools/sort/samtools/view_sam_to_bam/bowtie2/se_fasta_Abraham2017_hg19/edena/assembling_-d_20_-c_20_-minCoverage_5/edena/overlapping/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_H3K27ac_contigs.vcf.gz
    """
    input:
        txt="out/{filler}"
    output:
        txt="out/{tool}{extra}/{filler}"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="grep/"
    conda:
        "../envs/coreutils.yaml"
    shell:
        "grep {params.extra} {input.txt} > {output.txt}"
