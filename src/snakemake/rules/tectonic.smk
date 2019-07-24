rule tectonic_indir:
    """
    Created:
        2019-02-25 15:45:14
    Aim:
        Produce pdf from latex file. Tested as a replacement for texlive_pdflatex.
        This version compile in input directory then move the pdf to output directory. This is the simplest way to deal with tex file that depends on many other files, e.g. thesis template.
    Test:
        out/tectonic_indir/ln/updir/mw-gcthesis/src/tex/latexamu/these.pdf
    """
    input:
        tex="out/{filler}/{filename}.tex",
        dep=latex_smi_dep,
        #input_dep=latex_input_dependencies,
        #dep=latex_includegraphics_dependencies_test
    output:
        pdf="out/tectonic_indir/{filler}/{filename}.pdf"
    conda:
        "../envs/tectonic.yaml"
    log:
            "out/tectonic_indir/{filler}/{filename}.log"
    benchmark:
            "out/tectonic_indir/{filler}/{filename}.benchmark.tsv"
    threads:
        1
    shell:
        """
        (OUTDIR=`dirname {output.pdf}`
        INDIR=`dirname {input.tex}`
        for DIR in `ls -d */`;
        do
            ln -srf $DIR $INDIR/
            ln -srf $DIR $OUTDIR/
        done
        tectonic {input.tex}
        mv out/{wildcards.filler}/{wildcards.filename}.pdf {output.pdf}) &> {log}
        """

rule tectonic_softlink:
    """
    Created:
        2018-03-06 14:18:05
    Aim:
        Produce pdf from latex file. Tested as a replacement for texlive_pdflatex.
        This version was the first one but it produces longer log than "tectonic_indir"
        although it seems the output file from both rules are identical. 
        According to benchmarks, they have equivalent speed.
    Test:
        out/tectonic/ln/srf_from_inp/sbm/me/tex/invoice/facture-2018-10-001.pdf
        out/tectonic/ln/srf_from_src/tex/invoice/facture-2018-10-001.pdf
        out/tectonic/ln/srf_from_inp/me/tex/invoice/facture-2018-10-001.pdf
        out/tectonic/ln/updir/mw-gcthesis/src/tex/latexamu/these.pdf
    """
    input:
        tex="out/{filler}/{filename}.tex",
        dep=latex_smi_dep,
        #input_dep=latex_input_dependencies,
        #include_dep=latex_includegraphics_dependencies_test
    output:
        pdf="out/tectonic/{filler}/{filename}.pdf"
    log:
            "out/tectonic/{filler}/{filename}.log"
    benchmark:
            "out/tectonic/{filler}/{filename}.benchmark.tsv"
    conda:
        "../envs/tectonic.yaml"
    threads:
        1
    shell:
        """
        (OUTDIR=`dirname {output.pdf}`
        INDIR=`dirname {input.tex}`
        for DIR in `ls -d */`;
        do
            ln -srf $DIR $INDIR/
            ln -srf $DIR $OUTDIR/
        done
        tectonic --outdir $OUTDIR -p {input.tex}) &> {log}
        """

rule tectonic_softlink_no_dep:
    """
    Created:
        2018-03-06 14:18:05
    Aim:
        Produce pdf from latex file.
        Tested as a replacement for texlive_pdflatex.
        No dependencies are required here. 
        Use it in conjonction with \IfFileExists to get a compiled document with missing plots. 
        Useful when you need to preview the document and some plots are currently being produced.
    Test:
        "out/tectonic/beamer/nut_spike.pdf"
    """
    input:
        tex="out/{filler}.tex"
    output:
        pdf="out/tectonic_nodep/{filler}.pdf",
        tex="out/tectonic_nodep/{filler}.tex"
    threads:
        1
    conda:
        "../envs/tectonic.yaml"
    shell:
        """
        OUTDIR=`dirname {output.pdf}`

        ln -srf doc/ $OUTDIR/doc
        ln -srf inp/ $OUTDIR/inp
        ln -srf out/ $OUTDIR/out
        cp {input.tex} {output.tex}
        cd $OUTDIR
        tectonic `basename {output.tex}`

        #{input.tectonic} --outdir $OUTDIR {input.tex}
        """

rule tectonic_hardlink:
    """
    Created:
        2018-03-06 14:18:05
    Aim:
        Produce pdf from latex file. Tested as a replacement for texlive_pdflatex.
        This version of the rule ensure all plots in pdf are linked to original files in the same directory.
        Useful for plots you want to explore.
        Plots with long paths may not open on Windows due to FS path size limitations.
    Test:
        "out/tectonic_hardlink/beamer/sss/001_ChIP-seq_H3K27ac_Thymocytes/main.pdf"
    """
    input:
        tex="out/{filler}/{filename}.tex",
        input_dep=latex_input_dependencies,
        dep=latex_includegraphics_dependencies_test
    output:
        pdf="out/tectonic_hardlink/{filler}/{filename}.pdf",
        tex="out/tectonic_hardlink/{filler}/{filename}.tex"
    conda:
        "../envs/tectonic.yaml"
    threads:
        1
    shell:
        """
        OUTDIR=`dirname {output.pdf}`

        echo "OUTDIR: $OUTDIR"

        for DEP in {input}
        do
            echo "DEP: $DEP"
            LNDIR=`dirname $DEP`
            echo "LNDIR: $LNDIR"
            mkdir -p $OUTDIR/$LNDIR
            ln -fL $DEP $OUTDIR/$DEP
        done

        #cp {input.tex} {output.tex}
        cd $OUTDIR
        echo "moving to : `pwd`"
        tectonic `basename {output.tex}`
        """

rule tectonic_hardlink_no_dep:
    """
    Created:
        2018-10-19 22:06:46
    Aim:
        Produce pdf from latex file. Tested as a replacement for texlive_pdflatex.
        This version of the rule ensure all plots in pdf are linked to original files in the same directory.
        Useful for plots you want to explore.
        Plots with long paths may not open on Windows due to FS path size limitations.
        No dependencies are required here. 
        Use it in conjonction with \IfFileExists to get a compiled document with missing plots. 
        Useful when you need to preview the document and some plots are currently being produced.    Test:
        "out/tectonic_hardlink/beamer/sss/001_ChIP-seq_H3K27ac_Thymocytes/main.pdf"
    """
    input:
        tex="out/{filler}/{filename}.tex",
    output:
        pdf="out/tectonic_hardlink_no_dep/{filler}/{filename}.pdf",
        tex="out/tectonic_hardlink_no_dep/{filler}/{filename}.tex"
    conda:
        "../envs/tectonic.yaml"
    threads:
        1
    shell:
        """
        OUTDIR=`dirname {output.pdf}`
        echo "OUTDIR: $OUTDIR"
        cp {input.tex} {output.tex}
        cd $OUTDIR
        tectonic `basename {output.tex}`
        """


