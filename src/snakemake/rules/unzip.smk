rule unzip_p:
    """
    Created:
        2020-02-19 09:12:23
    Aim:
        Unzip only one zipped file.
    Test:
        out/unzip/p/wget/https/www.proteinatlas.org/download/proteinatlas.tsv
    """
    input:
        "out/{filler}.zip"
    output:
        "out/unzip/p/{filler}"
    shell:
        "unzip -p {input} > {output}"

rule unzip_d:
    """
    Created:
        2017-02-20 15:08:04
    Aim:
        Unzip zipped folder
    Test:
        "out/unzip/d/wget/ftp_ccb_jhu/pub/data/bowtie2_indexes/mm9/done"
        "out/unzip/d/wget/ftp/ftp.ccb.jhu.edu/pub/data/bowtie_indexes/hg19_1kgmaj_bt/done
    """
    input:
        "out/{filler}.zip"
    output:
        touch("out/unzip/d/{filler}/done")
    params:
        outdir="out/unzip/d/{filler}"
    shell:
        """
        unzip {input} -d {params.outdir}
        """

rule unzip_d_affymetrix_na36:
    input:
        "out/wget/http/www.affymetrix.com/Auth/analysis/downloads/na36/ivt/PrimeView.na36.annot.csv.zip"
    output:
        "out/unzip/d_affymetrix_na36/PrimeView.na36.annot.csv",
        "out/unzip/d_affymetrix_na36/3prime-IVT.AFFX_README.NetAffx-CSV-Files.txt"
    params:
        outdir="out/unzip/d_affymetrix_na36"
    shell:
        "unzip {input} -d {params.outdir}"

rule unzip_d_bowtie_index_hg19_1kgmaj:
    input:
        "out/wget/ftp/ftp.ccb.jhu.edu/pub/data/bowtie_indexes/hg19_1kgmaj_bt.zip"
    output:
         [x.strip() for x in open("../mw-lib/src/snakemake/lists/outputs_unzip_bowtie_index_hg19_1kgmaj.txt","r")]
    params:
        outdir="out/unzip/d_bowtie_index_hg19_1kgmaj"
    shell:
        "unzip {input} -d {params.outdir}"

rule unzip_d_bowtie_index_hg19:
    input:
        #"out/wget/ftp/ftp.ccb.jhu.edu/pub/data/bowtie_indexes/hg19_{1kgmaj}_bt.zip"
        "out/wget/ftp/ftp.ccb.jhu.edu/pub/data/bowtie_indexes/hg19.ebwt.zip"
    output:
         #[x.strip() for x in open("../mw-lib/src/snakemake/lists/outputs_unzip_bowtie_index_hg19_1kgmaj.txt","r")]
         expand("out/unzip/d_bowtie_index_hg19/hg19.{part}.ebwt", part=["1","2","3","4","rev.1","rev.2"])
    params:
        outdir="out/unzip/d_bowtie_index_hg19"
    shell:
        "unzip {input} -d {params.outdir}"

rule unzip_d_picard_tools_legacy:
    """
    Created:
        2017-02-20 15:08:04
    Test:
        out/unzip/d_picard_tools/picard-tools-1.99/picard-1.99.jar
        out/unzip/d_picard_tools/picard-tools-1.99/sam-1.99.jar
    """
    input:
        "out/wget/sourceforge/projects/picard/files/picard-tools/{version}/picard-tools-{version}.zip"
    output:
        picard="out/unzip/d_picard_tools/picard-tools-{version}/picard-{version}.jar",
        sam="out/unzip/d_picard_tools/picard-tools-{version}/sam-{version}.jar"
    wildcard_constraints:
        version="1.99" #versions below can be added here.
    params:
        outdir="out/unzip/d_picard_tools/"
    shell:
        """
        unzip {input} -d {params.outdir}
        """

rule unzip_d_miso_tools:
    """
    Created:
        2018-01-22 21:34:11
    Test:


        http://genes.mit.edu/burgelab/miso/annotations/miso_annotations_mm10_v1.zip
        http://genes.mit.edu/burgelab/miso/annotations/ver2/miso_annotations_mm10_v2.zip
        http://genes.mit.edu/burgelab/miso/annotations/ucsc_tables/mm10/ensGene.gff3

        extracting: out/unzip/d/wget/http/genes.mit.edu/burgelab/miso/annotations/miso_annotations_mm10_v1/mm10/log.txt
        extracting: out/unzip/d/wget/http/genes.mit.edu/burgelab/miso/annotations/miso_annotations_mm10_v1/mm10/TandemUTR.mm10.gff3
        extracting: out/unzip/d/wget/http/genes.mit.edu/burgelab/miso/annotations/miso_annotations_mm10_v1/mm10/SE.mm10.gff3
        extracting: out/unzip/d/wget/http/genes.mit.edu/burgelab/miso/annotations/miso_annotations_mm10_v1/mm10/RI.mm10.gff3
        extracting: out/unzip/d/wget/http/genes.mit.edu/burgelab/miso/annotations/miso_annotations_mm10_v1/mm10/MXE.mm10.gff3
        extracting: out/unzip/d/wget/http/genes.mit.edu/burgelab/miso/annotations/miso_annotations_mm10_v1/mm10/A3SS.mm10.gff3
        extracting: out/unzip/d/wget/http/genes.mit.edu/burgelab/miso/annotations/miso_annotations_mm10_v1/mm10/A5SS.mm10.gff3
        extracting: out/unzip/d/wget/http/genes.mit.edu/burgelab/miso/annotations/miso_annotations_mm10_v1/mm10/AFE.mm10.gff3
        extracting: out/unzip/d/wget/http/genes.mit.edu/burgelab/miso/annotations/miso_annotations_mm10_v1/mm10/ALE.mm10.gff3
        extracting: out/unzip/d/wget/http/genes.mit.edu/burgelab/miso/annotations/miso_annotations_mm10_v1/mm10/mm10_events_to_genes.zip

    """
    input:
        #"out/wget/sourceforge/projects/picard/files/picard-tools/{version}/picard-tools-{version}.zip"
    output:
        #picard="out/unzip/d_picard_tools/picard-tools-{version}/picard-{version}.jar",
        #sam="out/unzip/d_picard_tools/picard-tools-{version}/sam-{version}.jar"
    wildcard_constraints:
        version="1.99" #versions below can be added here.
    params:
        outdir="out/unzip/d_miso_tools/"
    shell:
        """
        unzip {input} -d {params.outdir}
        """

rule unzip_d_hmcan_diff:
    """
    Created:
        2018-02-25 18:16:36
    Test:

    """
    input:
        "out/wget/https/bitbucket.org/pyminer/hmcan-diff/get/3bd465fe4b76.zip"
    output:
        readme="out/unzip/d_hmcan_diff/pyminer-hmcan-diff-3bd465fe4b76/README.md",
        makefile="out/unzip/d_hmcan_diff/pyminer-hmcan-diff-3bd465fe4b76/src/Makefile"
    wildcard_constraints:
    params:
        outdir="out/unzip/d_hmcan_diff/"
    shell:
        """
        rm -rf out/unzip/d_hmcan_diff/
        unzip {input} -d {params.outdir}
        sed -i 's|-Wall|-Wall -I /gpfs/projects/spicuglia/mw/opt/miniconda/envs/hmcan_diff/include|g' {output.makefile}
        """




