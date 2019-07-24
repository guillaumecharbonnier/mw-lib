rule texlive_epstopdf:
    """
    Created:
        2017-06-08 15:51:33
    Aim:
        Convert eps file into pdf files.
    Warning:
        Current rule version works for test file but may fail for long path.
        epstopdf script is used inside this rule instead of directly by pdftolatex because it fails if provided too long paths for output pdf, typically around 114 characters, cf:
        http://computer-programming-forum.com/36-postscript/e29f2dadcc79639c.htm
    Test:
        out/texlive/epstopdf/deepTools/plotEnrichment_bed-mm10-test-srr-peaks_bam-mm10-test-srr.pdf
    """
    input:
        eps="out/{filler}/{sample}.eps",
    output:
        pdf="out/texlive/epstopdf/{filler}/{sample}.pdf"
    wildcard_constraints:
        sample="[^/]*"
    conda:
        "../envs/texlive_selected.yaml"
    threads:
        1
    shell:
        """
        epstopdf --outfile {output} {input}
        """
#
#        """
#        WDIR=`pwd`
#        OUTDIR=`dirname {output}`
#        cd $OUTDIR
#        epstopdf $WDIR/{input.eps} {wildcards.sample}.pdf
#        """

# Note: texlive pdftoeps does not exist, this may be the way to do this conversion:
#https://tex.stackexchange.com/questions/20883/how-to-convert-pdf-to-eps
rule texlive_pdftoeps:
    """
    Created:
        2019-02-15 11:40:38
    Not Working 
    Aim:
        Convert eps file into pdf files. epstopdf script is used inside this rule instead of directly by pdftolatex because it fails if provided too long paths for output pdf, typically around 114 characters, cf:
        http://computer-programming-forum.com/36-postscript/e29f2dadcc79639c.htm
    Test:
        out/texlive/pdftoeps/deepTools/plotEnrichment_bed-mm10-test-srr-peaks_bam-mm10-test-srr.eps
    """
    input:
        "out/{filler}.pdf",
    output:
        "out/texlive/pdftoeps/{filler}.eps"
    conda:
        "../envs/texlive_selected.yaml"
    threads:
        1
    shell:
        "pdftoeps {input} {output}"





