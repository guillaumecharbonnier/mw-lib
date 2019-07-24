"""
"""

#rule latex_report:
#    input:
#        pdflatex="opt/texlive/bin/x86_64-linux/pdflatex",
#        bibtex="opt/texlive/bin/x86_64-linux/bibtex",
#        tex="report/latex/report.tex",
#        code="report/latex/media/code/",
#        result="report/latex/media/result/",
#        request_figure_1="result/deepTools/plotHeatmap/groups/referencePoint/rep_pos/center/bamCoverage/mm9/for_article_4types/smt_max/1000bp/all.pdf"
#    output: 
#        pdf="report/latex/report.pdf"
#    threads: 1
#    shell:"""
#    WDIR=`pwd`
#    cd report/latex
#    $WDIR/{input.pdflatex} report.tex
#    $WDIR/{input.bibtex} report
#    # It is nice to compile again to get correct references.
#    $WDIR/{input.pdflatex} report.tex
#    
#    # Some plots are too big top include in the main document so every tex figure is compiled in a separate directory
#    mkdir --parents media/tex/fig/pdf
#    for fig in `ls media/tex/fig/*.tex`
#    do
#        $WDIR/{input.pdflatex} --output-directory media/tex/fig/pdf $fig
#    done
#    
#    # Tables
#    mkdir --parents media/tex/tab/pdf
#    for tab in `ls media/tex/tab/*.tex`
#    do
#        $WDIR/{input.pdflatex} --output-directory media/tex/tab/pdf $tab
#    done
#    
#    # Sections
#    mkdir --parents tex/sec/pdf
#    for sec in `ls tex/sec/*.tex`
#    do
#        $WDIR/{input.pdflatex} --output-directory tex/sec/pdf $sec
#    done
#    """

rule texlive_pdflatex:
    """
    Created:
        2016-07-03 13:07:49
    Modified:
        2017-01-16 10:25:12 - Files produced by latex are moved from 'report' to 'out/latex'. Input tex code is moved from 'report' to 'src/tex'.
        2017-08-01 10:58:20 - Removed {params.outdir} and replaced with OUTDIR which allows to call and compile correctly files in subdirectories of 'src/tex/', e.g. out/texlive/pdflatex/article/article_2.pdf.
    Aim:
        Produce pdf from latex file.
    Test:
        "out/texlive/pdflatex/update_2016_12_13.pdf"
    """
    input:
        tex="src/tex/{id}.tex",
        dep=latex_includegraphics_dependencies
    output:
        pdf="out/texlive/pdflatex/{id}.pdf"
    conda:
        "../envs/texlive_selected.yaml"
    threads:
        1
    shell:
        """
        OUTDIR=`dirname {output.pdf}`
        pdflatex --output-directory $OUTDIR {input.tex}
        """

rule get_rsync_incl_latex_includegraphics_dependencies:
    """
    Created:
        2017-07-26 11:08:58
    Aim:
        Get a file used by rsync to sync image dependencies for tex document in order to be able to compile from any host.
    Test:
        "src/rsync/incl_latex_includegraphics_dependencies/update_2016_12_13.txt"
        "src/rsync/incl_latex_includegraphics_dependencies/test_chromhmm.txt"

    Note:
        https://stackoverflow.com/questions/148451/how-to-use-sed-to-replace-only-the-first-occurrence-in-a-file
        0,/Apple/{s/Apple/Banana/}
        start at line 0, continue until you match 'Apple', execute the substitution in curly brackets. 
        
    """
    input:
        tex="src/tex/{id}.tex"
    params:
        dep=latex_includegraphics_dependencies
    output:
        txt="src/rsync/files-from/latex_includegraphics_dependencies/{id}.txt"
    threads:
        1
    shell:
        """
        echo '{params.dep}' | sed 's/ /\\n/g' > {output.txt}
        """

rule texlive_pdflatex_bibtex:
    """
    Created:
        2017-08-01 11:40:07
    Aim:
        Produce pdf from latex file using bibtex for citations.
    Test:
        "out/texlive/pdflatex_bibtex/article/test_report/report.pdf"
    """
    input:
        tex="src/tex/{id}.tex",
        dep=latex_includegraphics_dependencies
    output:
        pdf="out/texlive/pdflatex_bibtex/{id}.pdf"
    threads:
        1
    conda:
        "../envs/texlive_selected.yaml"
    shell:
        """
        INDIR=`dirname {input.tex}`
        OUTDIR=`dirname {output.pdf}`
        cp -r $INDIR/* $OUTDIR
        cd $OUTDIR
        ln -sf {WDIR}/out out
        ln -sf {WDIR}/doc doc
        pdflatex `basename {input.tex}`
        bibtex `basename {input.tex} | sed 's/.tex//'`
        pdflatex `basename {input.tex}`
        """

rule texlive_pdflatex_v2:
    """
    Created:
        2017-11-08 16:16:51
    Aim:
        Produce pdf from latex file using bibtex for citations.
    Test:
        "out/texlive/pdflatex_v2/invoice/2017-11-001.pdf"
    """
    input:
        tex="src/tex/{id}.tex",
        dep=latex_includegraphics_dependencies
    output:
        pdf="out/texlive/pdflatex_v2/{id}.pdf"
    conda:
        "../texlive_selected.yaml"
    threads:
        1
    shell:
        """
        INDIR=`dirname {input.tex}`
        OUTDIR=`dirname {output.pdf}`
        cp -r $INDIR/* $OUTDIR
        cd $OUTDIR
        ln -sf {WDIR}/out out
        ln -sf {WDIR}/doc doc
        pdflatex `basename {input.tex}`
        """


