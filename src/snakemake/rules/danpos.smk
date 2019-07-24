rule danpos_dtriple_extra:
    """
    Created:
        2019-01-22 17:11:46
    Doc:
          -m , --paired         set to 1 if the input data is mate-pair (paired-end) reads.  (default: 0)
    Test:
        out/danpos/dtriple/ln/alias/experiments/hg38_H3K27ac_LAL-T/CD34.wig -fn
    """
    input:
        bam="out/{filler}.bam"#,
        #python="opt/miniconda/envs/danpos/bin/python",
        #danpos="opt/danpos-2.2.2/danpos.py"
    output:
        wig =         "out/{tool}{extra}/{filler}/danpos.smooth.wig",
        ln_wig =      "out/{tool}{extra}/{filler}.wig",
        xls_positions="out/{tool}{extra}/{filler}/danpos.smooth.positions.xls",
        xls_peaks =   "out/{tool}{extra}/{filler}/danpos.smooth.peaks.xls",
        xls_regions = "out/{tool}{extra}/{filler}/danpos.smooth.regions.xls"
    params:
        extra = params_extra,
        ln_bam =      "out/{tool}{extra}/{filler}/danpos.bam",
        outdir =      "out/{tool}{extra}/{filler}/"
    wildcard_constraints:
        tool="danpos/dtriple"
    threads:
        1
    priority:
        2
    conda:
        "../envs/danpos.yaml"
    shell:
        """
        # Danpos include input file name in the output directory, which is ugly.
        # Link is used here so Danpos see clean input file.
        # Note: hardlink is bad with snakemake because when the hardlink is done on the bam there is a 'touch' on file, thus every rules using the bam have to be done again.
        WDIR=`pwd`
        ln -s --force $WDIR/{input.bam} $WDIR/{params.ln_bam}
        cd {params.outdir}

        INPUT_SAMPLE=`basename {params.ln_bam}`

        danpos.py dtriple {params.extra} --out . $INPUT_SAMPLE

        # Danpos put output files in './result/pooled/' when given '.' as '--out'.
        mv result/pooled/* ./
        rm -rf result
        ln $WDIR/{output.wig} $WDIR/{output.ln_wig}
        """

