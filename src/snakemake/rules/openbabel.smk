rule obabel_smi_to_ext:
    """
    Note:
        Created to convert smi to 2D structure. Usages could be broaden.
    Test:
        out/obabel/smi-to-svg/ln/updir/mw-gcthesis/inp/smi/5mC.svg
        out/obabel/smi-to-svg/ln/updir/mw-gcthesis/inp/smi/5hmC.svg
        out/obabel/smi-to-svg/ln/updir/mw-gcthesis/inp/smi/5fC.svg
        out/obabel/smi-to-svg/ln/updir/mw-gcthesis/inp/smi/5caC.svg
        out/obabel/smi-to-svg/ln/updir/mw-gcthesis/inp/smi/5hmU.svg
        out/obabel/smi-to-svg/ln/updir/mw-gcthesis/inp/smi/6mA.svg
    """
    input:
        smi="out/{filler}.{extin}"
    output:
        "out/{tool}{extin}-to-{extout}{extra}/{filler}.{extout}"
    wildcard_constraints:
        tool = "obabel/",
        extin = "smi|sdf",
        extout= "svg|png"
    params:
        extra = params_extra
    conda:
        "../envs/openbabel.yaml"
    shell:
        "babel {input} {output} {params.extra}"
