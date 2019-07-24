rule plantuml:
    """
    Created:
        2019-02-18 11:00:47
    Aim:
    Note:
        plantUML allows to produce svg, eps and pdf images but
        pdf fails and png aspect is way better than svg or eps
        as long as a decent dpi is set in the input file, e.g.:
        skinparam dpi 600

        Want to generate huge diagram? Check faq: http://plantuml.com/faq
    Test:
        out/plantuml/ln/updir/mw/src/plantuml/mapper.png
    """
    input:
        txt="out/{filler}.txt"
    output:
        txt="out/plantuml/{filler}.txt",
        png="out/plantuml/{filler}.png"
    conda:
        "../envs/plantuml.yaml"
    shell:
        "ln -srf {input.txt} {output.txt}; "
        "plantuml.sh -DPLANTUML_LIMIT_SIZE=65536 -Xmx2048m -tpng {output.txt}"

