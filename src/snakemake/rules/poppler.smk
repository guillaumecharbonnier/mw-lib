rule poppler_pdftoppm_png_singlefile:
    """
    Created:
        2017-03-16 11:51:53
    Aim:
        Tool to convert pdf to png. Especially to convert those vector plots so heavy they make slides and report lagging.
    Idea from:
        http://askubuntu.com/questions/50170/how-to-convert-pdf-to-image/50180#50180
    Test:
        out/poppler/pdftoppm_png_singlefile/deepTools/plotEnrichment_bed-mm10-test-srr-peaks_bam-mm10-test-srr.png
    """
    input:
        pdf="out/{filler}.pdf"
    output:
        png="out/poppler/pdftoppm_png_singlefile/{filler}.png"
    params:
        png_prefix="out/poppler/pdftoppm_png_singlefile/{filler}"
    conda:
        "../envs/poppler.yaml"
    shell:
        "pdftoppm -png -singlefile {input.pdf} {params.png_prefix}"
