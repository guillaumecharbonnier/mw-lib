"""
Interested in using some content of a multiple file archive in your workflow?
Here is how it can be done with snakemake:
1. Use 'rule tar_xvzf' to get an idea of the whole content of your archive.
2. Use 'ls' or 'find' to get the list of output files inside a text file:
    find $EXTRACTED_ARCHIVE -type f > src/snakemake/lists/archive_content.txt
3. Use this text file to define the output of a dedicated rule for your archive
4. Adjust paths in archive_content.txt to match the output directory of your dedicated rule.
sed 's|out/tar/xvzf/wget/ftp/ftp.ebi.ac.uk/pub/databases/blueprint/paper_data_sets/chromatin_states_carrillo_build37/BLUEPRINT_cell_lines|out/tar/xvzf_Carrillo2017_blueprint_cell_lines|' ../mw-lib/src/snakemake/lists/outputs_tar_xvzf_Carrillo2017_blueprint_cell_lines.txt -i

See for example 'rule tar_xvzf_Carrillo2017_roadmap'
"""

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
        out/tar/xvzf/wget/ftp/ftp.ebi.ac.uk/pub/databases/blueprint/paper_data_sets/chromatin_states_carrillo_build37/ROADMAP/done
        out/tar/xvzf/wget/ftp/ftp.ebi.ac.uk/pub/databases/blueprint/paper_data_sets/chromatin_states_carrillo_build37/BLUEPRINT_healthy/done
        out/tar/xvzf/wget/ftp/ftp.ebi.ac.uk/pub/databases/blueprint/paper_data_sets/chromatin_states_carrillo_build37/BLUEPRINT_cell_lines/done
        out/tar/xvzf/wget/ftp/ftp.ebi.ac.uk/pub/databases/blueprint/paper_data_sets/chromatin_states_carrillo_build37/BLUEPRINT_disease/done
        out/tar/xvzf/wget/ftp/ftp.ebi.ac.uk/pub/databases/blueprint/paper_data_sets/chromatin_states_carrillo_build37/ENCODE/done
        out/tar/xvzf/wget/https/cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A/done
    """
    input:
        "out/{filler}.tar.gz"
    output:
        touch("out/tar/xvzf/{filler}/done")
    params:
        outdir="out/tar/xvzf/{filler}"
    shell:
        "tar -xvzf {input} --directory {params.outdir}"

rule tar_xvzf_igenome:
    """
    Created:
        2017-05-23 11:18:05
    Aim:
        Extract gzipped tar archive. This rule is used for igenome archive where we expect index files in specific directories regardless of the species and assemblies.
        These archives also contain species- and assemblies- specific files (e.g. sequence of each chromosome) which can not be referenced by such generic rule.
        If you need to reference these specific files in your workflow, please see rule tar_xvzf_igenome_NCBI_GRCh38 below.
    Test:
        out/tar/xvzf_igenome/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome.1.bt2
        out/tar/xvzf_igenome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
        out/tar/xvzf_igenome/Homo_sapiens/NCBI/GRCh38/Sequence/Bowtie2Index/genome.1.bt2
    """
    input:
        # Testing for deprecation on 2018-01-10 12:18:41
        #"out/wget/ftp_igenome/{specie}/{source}/{index}/{specie}_{source}_{index}.tar.gz"
        "out/wget/ftp/igenome:G3nom3s4u@ussd-ftp.illumina.com/{specie}/{source}/{index}/{specie}_{source}_{index}.tar.gz"
    output:
        #expand("out/tar/xvzf_igenome/{{specie}}/{{source}}/{{index}}/{igenome_files}", igenome_files=["Sequence/Bowtie2Index/genome.1.bt2", "Sequence/WholeGenomeFasta/genome.fa", "Annotation/Genes/genes.gtf"])
        expand("out/tar/xvzf_igenome/{{specie}}/{{source}}/{{index}}/{igenome_files}", igenome_files=[x.strip() for x in open("../mw-lib/src/snakemake/lists/outputs_tar_xvzf_igenome.txt", "r")])
    params:
        outdir="out/tar/xvzf_igenome"
    wildcard_constraints:
        specie="[A-Za-z_]+",
        source="UCSC|Ensembl|NCBI"
    shell:
        """
        tar -xvzf {input} --directory {params.outdir}
        touch {output} # Without this touch the timestamp on the file makes it always older than the input tar archive, leading Snakemake to redo this rule every time.
        """

rule test_multiple_tar_xvzf_igenome:
    input:
        expand("out/tar/xvzf_igenome/{specie}/{source}/{index}/Sequence/Bowtie2Index/genome.1.bt2", zip, 
        specie=["Saccharomyces_cerevisiae","Drosophila_melanogaster","Mus_musculus","Mus_musculus","Homo_sapiens"],
        source=["UCSC"                    ,"UCSC"                   ,"Ensembl"     ,"UCSC"        ,"NCBI"],
        index =["sacCer3"                 ,"dm6"                    ,"GRCm38"      ,"mm9"         ,"GRCh38"])

rule tar_xvzf_igenome_NCBI_GRCh38:
    input: "out/wget/ftp/igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/NCBI/GRCh38/Homo_sapiens_NCBI_GRCh38.tar.gz"
    output: [x.strip() for x in open("../mw-lib/src/snakemake/lists/outputs_tar_xvzf_igenome_NCBI_GRCh38.txt", "r")]
    shell: "tar -xvzf {input} --directory out/tar/xvzf_igenome_NCBI_GRCh38"

rule tar_xvzf_Carrillo2017_roadmap:
    input: "out/wget/ftp/ftp.ebi.ac.uk/pub/databases/blueprint/paper_data_sets/chromatin_states_carrillo_build37/ROADMAP.tar.gz"
    output: [x.strip() for x in open("../mw-lib/src/snakemake/lists/outputs_tar_xvzf_Carrillo2017_roadmap.txt","r")]
    shell: "tar -xvzf {input} --directory out/tar/xvzf_Carrillo2017_roadmap"

rule tar_xvzf_Carrillo2017_encode:
    input: "out/wget/ftp/ftp.ebi.ac.uk/pub/databases/blueprint/paper_data_sets/chromatin_states_carrillo_build37/ENCODE.tar.gz"
    output: [x.strip() for x in open("../mw-lib/src/snakemake/lists/outputs_tar_xvzf_Carrillo2017_encode.txt","r")]
    shell: "tar -xvzf {input} --directory out/tar/xvzf_Carrillo2017_encode"

rule tar_xvzf_Carrillo2017_blueprint_healthy:
    input: "out/wget/ftp/ftp.ebi.ac.uk/pub/databases/blueprint/paper_data_sets/chromatin_states_carrillo_build37/BLUEPRINT_healthy.tar.gz"
    output: [x.strip() for x in open("../mw-lib/src/snakemake/lists/outputs_tar_xvzf_Carrillo2017_blueprint_healthy.txt","r")]
    shell: "tar -xvzf {input} --directory out/tar/xvzf_Carrillo2017_blueprint_healthy"

rule tar_xvzf_Carrillo2017_blueprint_cell_lines:
    input: "out/wget/ftp/ftp.ebi.ac.uk/pub/databases/blueprint/paper_data_sets/chromatin_states_carrillo_build37/BLUEPRINT_cell_lines.tar.gz"
    output: [x.strip() for x in open("../mw-lib/src/snakemake/lists/outputs_tar_xvzf_Carrillo2017_blueprint_cell_lines.txt","r")]
    shell: "tar -xvzf {input} --directory out/tar/xvzf_Carrillo2017_blueprint_cell_lines"

rule tar_xvzf_Carrillo2017_blueprint_disease:
    input:
        "out/wget/ftp/ftp.ebi.ac.uk/pub/databases/blueprint/paper_data_sets/chromatin_states_carrillo_build37/BLUEPRINT_disease.tar.gz"
    output:
        [x.strip() for x in open("../mw-lib/src/snakemake/lists/outputs_tar_xvzf_Carrillo2017_blueprint_disease","r")]
    shell:
        "tar -xvzf {input} --directory out/tar/xvzf_Carrillo2017_blueprint_disease"

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

rule tar_xvzf_genome_cellranger:
    input:
        "out/{filler}.tar.gz"
    output:
        touch("out/tar/xvzf_genome_cellranger/{filler}/done")
    params:
        outdir="out/tar/xvzf_genome_cellranger/{filler}"
    shell:
        "tar -xvzf {input} --directory `dirname {params.outdir}`"
