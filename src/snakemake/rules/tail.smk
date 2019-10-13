rule tail_extra:
    """
    Created:
        2019-09-09 16:23:49
    """
    input:
        "out/{filler}"
    output:
        "out/{tool}{extra}/{filler}"
    log:
        "out/{tool}{extra}/{filler}.log"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="tail/"
    shell:
        "tail {params.extra} {input} > {output} 2> {log}"
