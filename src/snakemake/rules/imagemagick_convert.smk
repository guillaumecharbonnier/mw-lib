rule imagemagick_convert_svg_to_png_extra:
    """
    Created:
        2017-05-10 10:52:12
    Aim:
        Convert svg to png
    Test:
        out/imagemagick/convert_svg-to-png_-trim/wget/apps_education_ucsb/w/images/3/3f/Todo_stickie.png
        out/imagemagick/convert_jpg-to-pdf/wget/https/marlin-prod.literatumonline.com/cms/attachment/3be4e7bd-9be8-4db5-8f38-ef77a6581270/gr1_lrg.pdf
        out/imagemagick/convert_svg-to-pdf/obabel/smi-to-svg/ln/updir/mw-gcthesis/inp/smi/6mA.pdf

    """
    input:
        image="out/{filler}.{extin}"
    output:
        image="out/{tool}_{extin}-to-{extout}{extra}/{filler}.{extout}"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "imagemagick/convert",
        # add other if necessary
        extin  = "png|jpg|pdf|svg",
        extout = "png|jpg|pdf|svg"
    conda:
        "../envs/imagemagick.yaml"
    shell:
        "convert {input.image} {params.extra} {output.image}"

rule imagemagick_convert_same_format_extra:
    """
    Created:
        2017-05-10 10:52:12
    Aim:
        Allow the use of convert transformation like "-trim" on an image file without converting its format.
    Test:
        out/imagemagick/convert-same-format_-trim/ChromHMM/LearnModel_test19_numstates-10_assembly-mm10/emissions_10.png
    """
    input:
        image="out/{filler}"
    output:
        image="out/{tool}{extra}/{filler}"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "imagemagick/convert-same-format"
    conda:
        "../envs/imagemagick.yaml"
    shell:
        "convert {input.image} {params.extra} {output.image}"
