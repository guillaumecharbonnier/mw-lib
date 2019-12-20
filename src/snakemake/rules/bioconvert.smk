rule bioconvert_inp2out:
    """
    Created:
        2018-11-21 16:17:43
    Aim:
        Convert xls files to tsv for easier view in terminal.
        Note: Never finished, I am trying with bioconvert instead.
    Test:
        out/bioconvert/xls2csv/wget/https/media.nature.com/original/nature-assets/ng/journal/v47/n10/extref/ng.3385-S4.csv
    """
    input:
        "out/{filler}.{inpext}"
    output:
        "out/{tool}{inpext}2{outext}{extra}/{filler}.{outext}"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "bioconvert/"
    conda:
        "../envs/bioconvert.yaml"
    shell:
        """
        bioconvert {wildcards.inpext}2{wildcards.outext} {params.extra} {input} {output}
        """

rule bioconvert_fastq2fasta_qual:
    """
    Created:
        2018-11-21 16:17:43
    Aim:
        Convert xls files to tsv for easier view in terminal.
        Note: Never finished, I am trying with bioconvert instead.
    Test:
        out/bioconvert/xls2csv/wget/https/media.nature.com/original/nature-assets/ng/journal/v47/n10/extref/ng.3385-S4.csv
    """
    input:
        "out/{filler}.fastq"
    output:
        temp("out/{tool}{extra}/{filler}.csfasta"),
        temp("out/{tool}{extra}/{filler}.qual")
    log:
        "out/{tool}{extra}/{filler}.log"
    benchmark:
        "out/{tool}{extra}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "bioconvert/fastq2fasta_qual"
    conda:
        "../envs/bioconvert.yaml"
    shell:
        "bioconvert fastq2fasta_qual {params.extra} {input} {output} &> {log}"

rule bioconvert_fasta_qual2fastq:
    """
    Created:
        2018-11-21 16:17:43
    Aim:
        Convert qual and fasta to fastq.
    Note:
        Current input suffixes correspond to output from SOLID but may be adapted for other usecases
    Test:
    """
    input:
        "out/{filler}.csfasta",
        "out/{filler}.QV.qual"
    output:
        temp("out/{tool}{extra}/{filler}.fastq")
    log:
        "out/{tool}{extra}/{filler}.log"
    benchmark:
        "out/{tool}{extra}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "bioconvert/fasta_qual2fastq"
    conda:
        "../envs/bioconvert.yaml"
    shell:
        "bioconvert fasta_qual2fastq {params.extra} {input} {output} &> {log}"

