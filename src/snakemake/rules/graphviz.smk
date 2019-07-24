rule graphviz_dot:
    """
    Created:
        2019-02-08 16:23:36
    Aim:
        out/graphviz/dot/snakemake/stdout_--rulegraph_test.pdf
    """
    input:
        "out/{filler}.dot"
    output:
        "out/{tool}{extra}/{filler}.{ext}"
    params:
        extra = params_extra
    conda:
        "../envs/graphviz.yaml"
    wildcard_constraints:
        tool="graphviz/dot"
    shell:
        "dot {input} -T{wildcards.ext} -o {output} {params.extra}"
