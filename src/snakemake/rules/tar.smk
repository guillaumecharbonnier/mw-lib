rule tar_xvzf:
    """
    Created:
        2017-02-13 14:02:42
    Aim:
        Extract gzipped tar archive. 
    Test:
        out/tar/xvzf/input/atlas-latest-data
        out/tar/xvzf/wget/ftp_ncbi/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index

        out/tar/xvzf/wget/ftp_igenome/Saccharomyces_cerevisiae/UCSC/sacCer3/Saccharomyces_cerevisiae_UCSC_sacCer3
        out/tar/xvzf/wget/ftp_igenome/Drosophila_melanogaster/UCSC/dm6/Drosophila_melanogaster_UCSC_dm6
        out/tar/xvzf/wget/http/www.regulatory-genomics.org/wp-content/uploads/2015/07/THOR-tools
        out/tar/xvzf/wget/http/www.cbrc.kaust.edu.sa/hmcan/hmcan-diff_example
    """
    input:
        "out/{filler}.tar.gz"
    output:
        touch("out/tar/xvzf/{filler}/done")
    params:
        outdir="out/tar/xvzf/{filler}"
    shell:
        """
        tar -xvzf {input} --directory {params.outdir}
        """

rule tar_xvzf_igenome:
    """
    Created:
        2017-05-23 11:18:05
    Aim:
        Extract gzipped tar archive. This rule is used for igenome archive where we expect bowtie index in specific directories.
    Test:
        out/tar/xvzf_igenome/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index
        out/tar/xvzf_igenome/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome.1.bt2
        input:
            out/wget/ftp/igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/UCSC/hg38/Homo_sapiens_UCSC_hg38.tar.gz
        output:
            out/tar/xvzf_igenome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
    Note:
        More output files could be explicited here if needed.
    """
    input:
        # Testing for deprecation on 2018-01-10 12:18:41
        #"out/wget/ftp_igenome/{specie}/{source}/{index}/{specie}_{source}_{index}.tar.gz"
        "out/wget/ftp/igenome:G3nom3s4u@ussd-ftp.illumina.com/{specie}/{source}/{index}/{specie}_{source}_{index}.tar.gz"
    output:
        "out/tar/xvzf_igenome/{specie}/{source}/{index}/Sequence/Bowtie2Index/genome.1.bt2",
        "out/tar/xvzf_igenome/{specie}/{source}/{index}/Sequence/WholeGenomeFasta/genome.fa",
        "out/tar/xvzf_igenome/{specie}/{source}/{index}/Annotation/Genes/genes.gtf"
    params:
        outdir="out/tar/xvzf_igenome"
    wildcard_constraints:
        specie="[A-Za-z_]+",
        source="UCSC|Ensembl|NCBI"
    shell:
        """
        mkdir -p {output}
        tar -xvzf {input} --directory {params.outdir}
        touch {output} # Without this touch the timestamp on the file makes it always older than the input tar archive, leading Snakemake to redo this rule every time.
        """

rule tar_xvzf_danpos:
    """

    """
    input:
        "out/wget/http/lilab.research.bcm.edu/dldcc-web/lilab/kaifuc/danpos/release/danpos-2.2.2.tgz"
    output:
        scripts=expand("opt/danpos-2.2.2/danpos.py", script=["danpos","wig","wigs","wig"])
    params:
        outdir="opt"
    shell:
        """
        mkdir -p {params.outdir}
        tar -xvzf {input} --directory {params.outdir}
        """

rule tar_xvf_data_brdt:
    """
    Todo:
        Generalize the rule for all .tar files. Or just let it deprecated and use directly tar with -z argument...
    Test:
        "out/gunzip/to-stdout/tar/xvf_data_brdt/GSM984200_4_EM2_R_BRDT-Paired.F3.no_4-W200-G600-islands-summary-FDR001.bed"
    """
    input:
        "out/wget/ftp_ncbi/geo/series/GSE39nnn/GSE39910/suppl/GSE39910_RAW.tar"
    output:
        expand("out/tar/xvf_data_brdt/{id}",
            id=[
                "GSM2101805_7_P_BD1_F3-F5-BC-Paired2-W200-G600-islands-summary-FDR001.bed.gz",
                "GSM2101806_8_R_BD1_F3-F5-BC-Paired2-W200-G600-islands-summary-FDR001.bed.gz",
                "GSM984199_3_EM2_P_BRDT-Paired.F3.no_4-W200-G600-islands-summary-FDR001.bed.gz",
                "GSM984200_4_EM2_R_BRDT-Paired.F3.no_4-W200-G600-islands-summary-FDR001.bed.gz"
                ])
    params:
        outdir="out/tar/xvf_data_brdt"
    shell:
        """
        WDIR=`pwd`
        mkdir -p {params.outdir}
        cd {params.outdir}
        tar xvf $WDIR/{input}
        """

rule tar_xvzf_hmcan_diff_test_example:
    """
    Created:
        2018-02-26 13:56:44
    Aim:
        Extract gzipped tar archive. 
    Test:
        out/tar/xvzf_hmcan_diff_test_example/hmcan-diff_example/control_C1.sam
    Note:
        hmcan-diff_example/
        hmcan-diff_example/control_C1.sam
        hmcan-diff_example/chip_C2_rep_3.sam
        hmcan-diff_example/chip_C1_rep_3.sam
        hmcan-diff_example/C1_files.txt
        hmcan-diff_example/C2_files.txt
        hmcan-diff_example/chip_C1_rep_2.sam
        hmcan-diff_example/C2_control.txt
        hmcan-diff_example/chip_C2_rep_2.sam
        hmcan-diff_example/C1_control.txt
        hmcan-diff_example/chip_C1_rep_1.sam
        hmcan-diff_example/chip_C2_rep_1.sam
        hmcan-diff_example/reference/
        hmcan-diff_example/reference/chr1.fa
        hmcan-diff_example/data/
        hmcan-diff_example/data/hg19-blacklist.bed
        hmcan-diff_example/data/GC_profile_100KbWindow_Mapp76_hg19.cnp
        hmcan-diff_example/control_C2.sam
    """
    input:
        "out/wget/http/www.cbrc.kaust.edu.sa/hmcan/hmcan-diff_example.tar.gz"
    output:
        expand("out/tar/xvzf_hmcan_diff_test_example/hmcan-diff_example/{filler}", filler=TAR_CONTENT_HMCAN_DIFF_TEST_EXAMPLE)
    params:
        outdir="out/tar/xvzf_hmcan_diff_test_example/"
    shell:
        """
        mkdir -p {output}
        tar -xvzf {input} --directory {params.outdir}
        """


