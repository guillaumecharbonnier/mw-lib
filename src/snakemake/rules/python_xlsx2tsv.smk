rule python_xlsx2tsv:
    """
    Created:
        2018-01-25 11:08:57
    Aim:
        Convert xlsx files to tsv for easier view in terminal, and use as input for Salva workflow
    Test:
        out/python/xlsx2tsv/wget/salva_runs.tsv
        out/python/xls2tsv/wget/https/media.nature.com/original/nature-assets/ng/journal/v47/n10/extref/ng.3385-S4.tsv


        src/python/xlsx2tsv.py
    """
    input:
        xlsx2tsv="../mw-lib/src/python/xlsx2tsv.py",
        xlsx="out/{filler}.xlsx"
    output:
        tsv="out/python/xlsx2tsv/{filler}.tsv"
    conda:
        "../envs/openpyxl.yaml"
    shell:
        "{input.xlsx2tsv} {input.xlsx} > {output.tsv}"

rule python_xls2tsv:
    """
    Created:
        2018-11-21 16:08:10
    Aim:
        Convert xls files to tsv for easier view in terminal.
        Note: Never finished, I am trying with bioconvert instead.
    Test:
        out/python/xls2tsv/wget/https/media.nature.com/original/nature-assets/ng/journal/v47/n10/extref/ng.3385-S4.tsv
    """
    input:
        xls2tsv="../mw-lib/src/python/xls2tsv.py",
        xls="out/{filler}.xls"
    output:
        tsv="out/python/xls2tsv/{filler}.tsv"
    conda:
        "../envs/xlrd.yaml"
    shell:
        "{input.xls2tsv} {input.xls} > {output.tsv}"

