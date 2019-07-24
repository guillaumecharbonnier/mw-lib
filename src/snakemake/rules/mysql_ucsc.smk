rule mysql_extra:
    """
    Created:
        2019-02-15 19:11:13
    Aim:
        Use mysql to download various tabulated files
    Test:
        out/mysql/ucsc-cpgIslandExt-mm10.bed
    """
    output:
        "out/{tool}{extra}"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="mysql/"
    conda:
        "../envs/mysql.yaml"
    shell:
        "mysql {params} > {output}"

rule mysql_ucsc_get_cpgIslandExt:
    """
    Created:
        2016-12-02 11h20
    Test:
        out/mysql/ucsc/mm10/cpgIslandExt.bed
    """
    output:
        bed="out/mysql/ucsc/{index}/cpgIslandExt.bed"
    conda:
        "../envs/mysql.yaml"
    shell:
        """
        mysql --user=genome -N \
        --host=genome-mysql.cse.ucsc.edu \
        -A -D {wildcards.index} \
        -e "select chrom,chromStart,chromEnd,name,length,cpgNum,gcNum,perCpg,perGc,obsExp from cpgIslandExt" > {output.bed}
        """

rule mysql_ucsc_get_chromInfo:
    """
    Tips: 
        One can use the UCSC Genome Browser's MySQL database to extract
        chromosome sizes. For example, H. sapiens:

        mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
        "select chrom, size from hg19.chromInfo"  > hg19.genome
    Test:
        out/mysql/ucsc/mm10/chromInfo_all.tsv
    """
    output:
        tsv_all  = "out/mysql/ucsc/{index}/chromInfo_all.tsv",
        tsv_main = "out/mysql/ucsc/{index}/chromInfo_main.tsv"
    conda:
        "../envs/mysql.yaml"
    shell:
        """
        mysql --user=genome -N \
        --host=genome-mysql.cse.ucsc.edu \
        -A -D {wildcards.index} \
        -e "select chrom, size from chromInfo" > {output.tsv_all}
    
        sed '/random/d' {output.tsv_all} > {output.tsv_main}
        """

rule mysql_ucsc_get_ensGene_D:
    """
    Created:
        2016-04-04 18h45
    Modified:
        2017-03-21 16:14:55
    Test:
        out/mysql/ucsc_get_ensGene_D-mm9/ensGene.bed
    """   
    output:
        ensgene="out/mysql/ucsc_get_ensGene_D-{index}/ensGene.bed"
    conda:
        "../envs/mysql.yaml"
    shell:
        """
        mysql --user=genome -N \
        --host=genome-mysql.cse.ucsc.edu \
        -A -D {wildcards.index} \
        -e "select chrom,txStart,txEnd,name,strand,name2 from ensGene" > {output.ensgene}
        """

