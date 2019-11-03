rule datamash_transpose:
    """
    out/datamash/transpose/bedtools/unionbedg_-header/bed-hg19-segments-thymopoiesis-tall-samples.bg

    """
    input:
        "out/{filler}"
    output:
        "out/{tool}{extra}/{filler}"
    wildcard_constraints:
        tool="datamash/"
    params:
        extra = params_extra
    conda:
        "../envs/datamash.yaml"
    shell:
        "datamash {params.extra} < {input} > {output}"

