rule sample_extra:
    """
    Created:
        2018-11-03 18:09:00
    Aim:
        Sample lines with better memory management than "shuf" hence more appropriate for huge bed files like small structures from MNase-Seq.
    Test:
        out/sample/_-k_10000000/
    """
    input:
        "out/{filler}"
    output:
        "out/{tool}{extra}/{filler}"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "sample/"
    conda:
        "../envs/sample.yaml"
    shell:
        "sample {params.extra} {input} > {output}"
