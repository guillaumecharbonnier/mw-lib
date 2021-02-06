rule danpos_wiq:
    """
    Created:
        2018-01-12 17:15:32
    Aim:
        Testing Danpos wig to wiq conversion.
    Test:
        dtriple_width-10_distance-60_edge-1_paired-1/

        wig =         "out/danpos/dtriple_width-{width}_distance-{distance}_edge-{edge}_paired-{paired}/{filler}/danpos.smooth.wig",
        out/danpos/dtriple_width-10_distance-100_edge-1_paired-0/inp/bam/hg38/H3K27ac/thymus/Blueprint/Th101_LC/danpos.smooth.wig
        out/danpos/wig_/danpos/dtriple_width-10_distance-100_edge-1_paired-0/inp/bam/hg38/H3K27ac/thymus/Blueprint/Th101_LC_QnormVS_Th91_LC.wig
    out/danpos/wiq_hg38/danpos/dtriple_width-10_distance-100_edge-1_paired-0/inp/bam/hg38/H3K27ac/thymus/Blueprint/Th101_LC_QnormVS_Th118_LC.wiq
    Note:
        Th125_SP8 seems to be a good reference based on ACTB and GAPDH signal to noise ratio in IGV.
        out/danpos/wiq_hg38/danpos/dtriple_width-10_distance-100_edge-1_paired-0/inp/bam/hg38/H3K27ac/thymus/Blueprint/Th101_LC_qnorVS_Th118_LC/Th101_LC.wig
    """
    input:
        wig_reference = "out/{filler}/{reference}.wig",
        wig_to_wiq = "out/{filler}/{wig_to_wiq}.wig",
        chrominfo = lambda wildcards: eval(mwconf['ids'][wildcards.chrominfo_id])
    output:
        #wig = "out/danpos/wiq_{chrominfo_id}/{filler}/{wig_to_wiq}_qnorVS_{reference}/{wig_to_wiq}.wig"
        wig = "out/danpos/wiq_{chrominfo_id}/{filler}/{wig_to_wiq}_qnorVS_{reference}.wig"
    log:
        "out/danpos/wiq_{chrominfo_id}/{filler}/{wig_to_wiq}_qnorVS_{reference}.log"
    benchmark:
        "out/danpos/wiq_{chrominfo_id}/{filler}/{wig_to_wiq}_qnorVS_{reference}.benchmark.tsv"
    params:
        outdir = "out/danpos/wiq_{chrominfo_id}/{filler}/{wig_to_wiq}_qnorVS_{reference}"
    conda:
        "../envs/danpos.yaml"
    wildcard_constraints:
        wig_to_wiq = "[a-zA-Z0-9_-]+",
        reference = "[a-zA-Z0-9_-]+"
    threads:
        1
    priority:
        2
    shell:
        """
        # Danpos include input file name in the output directory, which is ugly.
        # Link is used here so Danpos see clean input file.
        # Note: hardlink is bad with snakemake because when the hardlink is done on the bam there is a 'touch' on file, thus every rules using the bam have to be done again.
        ( WDIR=`pwd`
        mkdir -p {params.outdir}
        ln -srf {input.wig_to_wiq} {params.outdir}/wig_to_wiq.wig
        ln -srf {input.wig_reference} {params.outdir}/reference.wig
        ln -srf {input.chrominfo} {params.outdir}/chrominfo
        cd {params.outdir}
        danpos.py wiq chrominfo wig_to_wiq.wig --reference reference.wig

        ln -f wiq_result/wig_to_wiq.qnor.wig $WDIR/{output.wig}

        #TODO: remove this comment when the rule is working:
        rm -rf $WDIR/{params.outdir} ) &> {log}
        """
