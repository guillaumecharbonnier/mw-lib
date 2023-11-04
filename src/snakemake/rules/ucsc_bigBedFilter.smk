rule ucsc_bigBedFilter:
    """
    Created:
        2023-07-01
    Aim:
        I got the JASPAR2022.bb file from UCSC, but it does not load in igv-webapp (111GB).
        Here I try to filter it to only include the motifs with the highest score
    Test:
        out/ucsc/bigBedFilter_-minScore_400/JASPAR2022.bb
    """
    input:
        bb = "out/{filler}.bb"
    output:
        bb="out/{tool}{extra}/{filler}.bb"
    log:
        "out/{tool}{extra}/{filler}.log"
    benchmark:
        "out/{tool}{extra}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="ucsc/bigBedFilter"
    conda:
        "../envs/ucsc.yaml"
    shell:
        "bigBedFilter  {input.bb} {params.extra} {output.bb} &> {log}"
