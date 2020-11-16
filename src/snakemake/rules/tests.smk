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
        star_pe_extra="out/star/pe_fastq.gz_to_bam_staridx-GRCh38-ensembl_gtf-GRCh38-ensembl/sickle/pe_-t_sanger_-q_30/sra-tools/fastq-dump_pe/SRR2040277.bam
",
