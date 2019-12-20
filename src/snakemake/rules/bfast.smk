rule bfast_solid2fastq:
    """
    Created:
        2019-12-10 12:08:44
    Aim:
    Note:
    Test:
    """
    input:
        csfasta="out/{filler}.csfasta",
        qual="out/{filler}.QV.qual"
    output:
        "out/{tool}{extra}/{filler}.fastq"
    log:
        "out/{tool}{extra}/{filler}.log"
    benchmark:
        "out/{tool}{extra}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "bfast/solid2fastq"
    conda:
        "../envs/bfast.yaml"
    shell:
        "solid2fastq -o out/{wildcards.tool}{wildcards.extra}/{wildcards.filler} {params.extra} {input} &> {log}"


