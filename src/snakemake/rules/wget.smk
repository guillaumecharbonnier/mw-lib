rule wget_protocol:
    """
    Created:
        2017-09-07 15:12:23
    Modified:
        2018-01-10 11:55:33 - Extended to https and ftp.
    Aim:
        Download any files available from internet. If the needed URL contains special characters or is too long to be kept as a convenient filepath, you can use wget_extra instead.
   Test:
        out/wget/ftp/ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/SRA_Accessions.tab
        out/wget/http/pmgenomics.ca/hoffmanlab/proj/segway/2011/test.genomedata
        out/wget/https/www.tau.ac.il/~elieis/HKG/HK_genes.txt
        out/wget/ftp/ftp.ensembl.org/pub/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh38.87.chr.gtf.gz
        out/wget/ftp/igenome:G3nom3s4u@ussd-ftp.illumina.com/Drosophila_melanogaster/UCSC/dm6/Drosophila_melanogaster_UCSC_dm6.tar.gz
        out/wget/ftp/igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/UCSC/hg38/Homo_sapiens_UCSC_hg38.tar.gz
    """
    output:
        "out/wget/{protocol}/{path}"
    log:
        "out/wget/{protocol}/{path}.log"
    benchmark:
        "out/wget/{protocol}/{path}.benchmark.tsv"
    wildcard_constraints:
        protocol="http|https|ftp"
    shell:
        "wget --output-file={log} --output-document={output} {wildcards.protocol}://{wildcards.path}"

rule wget_extra:
    """
    Created:
        2019-02-05 23:41:55
    Aim:
        wget_protocol is an extremely convenient rule. However, some URL may contain special characters. In this case wget_extra allows to declare the URL in src/snakemake/tables/extraids.tsv and use an ID in the filepath instead.
    Test:
        out/wget/hg38-GENCODE-knownGene.bed
    """
    output:
        "out/{tool}{extra}"
    log:
        "out/{tool}{extra}.log"
    benchmark:
        "out/{tool}{extra}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="wget/"
    shell:
        "wget --output-file={log} --output-document={output} {params.extra}"

