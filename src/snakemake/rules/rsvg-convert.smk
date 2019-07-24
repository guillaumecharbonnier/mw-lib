rule rsvg_convert_extra:
    """
    Created:
        2018-11-12 13:09:54
    Note:
        Tested but produce empty pdf:
        out/rsvg-convert/svgtopdf/ChromHMM/LearnModel_test19_numstates-10_assembly-mm10/emissions_10.pdf
        out/rsvg-convert/svgtopdf/plantuml/ln/updir/mw/src/plantuml/test.pdf

        out/rsvg-convert/svgtopdf/wget/https/upload.wikimedia.org/wikipedia/commons/4/48/Animal_cell_structure_en.pdf
    """
    input:
        svg="out/{filler}.svg"
    output:
        pdf="out/{tool}{extra}/{filler}.pdf"
    wildcard_constraints:
        tool = "rsvg-convert/svgtopdf"
    params:
        extra = params_extra
    conda:
        "../envs/r_rsvg.yaml"
    shell:
        "rsvg-convert {params.extra} -f pdf -o {output.pdf} {input.svg}"
