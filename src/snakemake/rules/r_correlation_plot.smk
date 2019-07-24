rule correlation_plot_bw2matrix:
    """
    Modified:
        2016-04-26 9h56 - Adapted from part1.
    Aim:
        This rule take a matrix produced by bw2Matrix and create correlation plots.
    Usage:
        expand("result/corrPlot/{index}/{id}/{correlationType}/done", index="mm9", id="all_psk_subsamples", correlationType=["spearman","pearson"])
    """
    input:
        rscript="opt/anaconda3/bin/Rscript",\
        code="code/r/corrPlot.R",\
        tsv="out/bw2Matrix/{signal}/{index}/{id}.tsv"
    output:
        pdf="result/corrPlot/{signal}/{index}/{id}/{correlationType}/CorrPlot_{correlationType}.pdf"
    threads: 1
    shell:"""
    {input.rscript} {input.code} -i {input.tsv} -o result/corrPlot/{wildcards.signal}/{wildcards.index}/{wildcards.id}/{wildcards.correlationType} -c {wildcards.correlationType}
    """

rule correlation_plot_gtftk_coverage_matrix:
    """
    Created:
        2017-03-10 10:16:14
    Aim:
        This rule take a matrix produced by gtftk coverage with --matrix argument and create correlation plots.
    Test:
        "out/r/correlation_plot_{correlationType}/{filler}.pdf"
    """
    input:
        rscript="opt/miniconda/envs/bin/Rscript",
        code="code/r/script/corrPlot.R",
        txt="out/{filler}.txt"
    output:
        pdf="out/r/correlation_plot_{correlationType}/{filler}.pdf"
    params:
        outdir="out/r/correlation_plot_{correlationType}/{filler}"
    wildcard_constraints:
        correlationType="pearson|spearman|kendall"
    threads:
        1
    shell:
        """
        {input.rscript} {input.code} \
            -i {input.txt} \
            -o {params.outdir} \
            -c {wildcards.correlationType}
        """

