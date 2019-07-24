rule inkscape_svgtopdf_extra:
    """
    Created:
        2018-11-12 12:59:57
    Note:
        Inkscape currently unavailable on conda Linux Channel. Trying cairosvg as a replacement.
    """
    input:
        svg="out/{filler}.svg"
    output:
        pdf="out/{tool}{extra}/{filler}.pdf"
    wildcard_constraints:
        tool = "inkscape/svgtopdf"
    params:
        extra = params_extra
    conda:
        "../envs/inkscape.yaml"
    shell:
        "inkscape {params.extra} {input.svg} --export-pdf={output.pdf}"


