"""
Extract tests string
sed -n '/Test:/,/""/p' *.rules | sed -n '/Test:/!p' | sed -n '/""/!p'
"""

rule light_unit_tests:
    input:
        wget="out/wget/ftp/ftp.ensembl.org/robots.txt",

rule heavy_unit_tests:
    """
    Aim:
        Tests for rules with huge resources need.
        Only DAG generation is tested in CI.
    """
    input:
        sra_tools_fastq_dump_se_extra="out/sra-tools/fastq-dump_se/SRR1202037.fastq",

rule integration_tests:
    input:
        star_pe_extra="out/star/pe_fastq.gz_to_bam_staridx-GRCh38-ensembl_gtf-GRCh38-ensembl/sickle/pe_-t_sanger_-q_30/sra-tools/fastq-dump_pe/SRR2040277.bam",

rule cache_tests:
    """
    Aim:
        The [experimental cache directive](https://snakemake.readthedocs.io/en/stable/executing/caching.html) can overwrite different jobs of the same rule because it does not include the output path in the hash. Here I test that outputs for rules using cache are correctly written to avoid this confusion.
        These 4 files should have a different md5sum
    """
    input:
        "out/wget/ftp/ftp.ensembl.org/pub/release-102/removed_genomes.txt",
        "out/wget/ftp/ftp.ensembl.org/pub/release-102/renamed_genomes.txt",
        "out/wget/hg38-knownGene.bed",
        "out/wget/GSE231300-ENCFF826CKB-conservative-IDR-thresholded-peaks-GRCh38.bed.gz"
    shell:
        "md5sum {input}"

