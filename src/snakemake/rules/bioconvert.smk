rule bioconvert_inp2out:
    """
    Created:
        2018-11-21 16:17:43
    Aim:
        Convert xls files to tsv for easier view in terminal.
        Note: Never finished, I am trying with bioconvert instead.
    Test:
        out/bioconvert/xls2csv/wget/https/media.nature.com/original/nature-assets/ng/journal/v47/n10/extref/ng.3385-S4.csv
    """
    input:
        "out/{filler}.{inpext}"
    output:
        "out/{tool}{inpext}2{outext}{extra}/{filler}.{outext}"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "bioconvert/"
    conda:
        "../envs/bioconvert.yaml"
    shell:
        """
        bioconvert {wildcards.inpext}2{wildcards.outext} {params.extra} {input} {output}
        """


