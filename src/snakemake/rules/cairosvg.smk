rule cairosvg_svgtopdf_extra:
    """
    Created:
        2018-11-12 13:09:54
    Note:
        Tested but produce empty pdf:
        out/cairosvg/svgtopdf/ChromHMM/LearnModel_test19_numstates-10_assembly-mm10/emissions_10.pdf
    """
    input:
        svg="out/{filler}.svg"
    output:
        pdf="out/{tool}{extra}/{filler}.pdf"
    wildcard_constraints:
        tool = "cairosvg/svgtopdf"
    params:
        extra = params_extra
    conda:
        "../envs/cairosvg.yaml"
    shell:
        "cairosvg {params.extra} {input.svg} -o {output.pdf}"


