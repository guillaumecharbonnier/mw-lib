rule gs_compress_extra:
    """
    Created:
        2018-11-07 05:54:05
    """
    input:
        pdf="out/{filler}.pdf"
    output:
        pdf="out/{tool}{extra}/{filler}.pdf"
    wildcard_constraints:
        tool = "gs/compress"
    params:
        extra = params_extra
    shell:
        "gs {params.extra} -sOutputFile={output.pdf} {input.pdf}"

rule gs_sDEVICE_dCompatibilityLevel_dPDFSETTINGS_legacy:
    """
    Created:
        2018-08-20 18:09:05
    Aim:
        This is legacy gs_compress_extra.
        Various operations using gs, especially for pdf manipulation.
    Examples:
        Compress using dPDFSETTINGS="/screen":
        Set pdf version to 1.5 because when includepdf only work for pdf versions up to 1.5:
            out/gs/sDEVICE-pdfwrite_dCompatibilityLevel-1.5_dPDFSETTINGS-prepress/ln/alias/mendeley_path_to_bibtex_key/Goudarzi2016.pdf
        out/gs/sDEVICE-pdfwrite_dCompatibilityLevel-1.5_dPDFSETTINGS-screen/doc/inkscape/id/carte_identite_guillaume_charbonnier.pdf
    Notes:
        -dPDFSETTINGS=/screen   (screen-view-only quality, 72 dpi images)
        -dPDFSETTINGS=/ebook    (low quality, 150 dpi images)
        -dPDFSETTINGS=/printer  (high quality, 300 dpi images)
        -dPDFSETTINGS=/prepress (high quality, color preserving, 300 dpi imgs)
        -dPDFSETTINGS=/default  (almost identical to /screen)
    Test:
        out/gs/sDEVICE-pdfwrite_dCompatibilityLevel-1.5_dPDFSETTINGS-prepress/tectonic_hardlink/ln/srf_from_src/tex/beamer/sss/002_Integrative_Omics_Thymocytes/supplementary_figures.pdf
    """
    input:
        pdf="out/{filler}.pdf"
    output:
        pdf="out/gs/sDEVICE-{sDEVICE}_dCompatibilityLevel-{dCompatibilityLevel}_dPDFSETTINGS-{dPDFSETTINGS}/{filler}.pdf"
    wildcard_constraints:
        sDEVICE="pdfwrite",
        dCompatibilityLevel="1.4|1.5",
        dPDFSETTINGS="screen|ebook|printer|prepress|default"
    shell:
        """
        gs \
            -dNOPAUSE \
            -dBATCH \
            -sDEVICE={wildcards.sDEVICE} \
            -dCompatibilityLevel={wildcards.dCompatibilityLevel} \
            -dPDFSETTINGS="/{wildcards.dPDFSETTINGS}" \
            -sOutputFile={output.pdf} \
            {input.pdf}
        """

