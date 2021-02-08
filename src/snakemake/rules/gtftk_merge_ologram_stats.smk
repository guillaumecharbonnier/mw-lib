rule gtftk_merge_ologram_stats:
    """
    Created:
        2019-04-02 16:19:18
    Aim:
        OverLap Of Genomic Regions Analysis using Monte Carlo to annotate peaks.
    Doc:
        https://dputhier.github.io/pygtftk/annotation.html#ologram
    Test:
        out/gtftk/ologram_tfbs_chrominfo-hg19_gtf-hg19-tfbs/bedtools/merge_-d_5/sort/_-k1,1_-k2,2n/crossmap/chain-hg38-to-hg19/r/order_cd34_ec_common_atac_peaks_according_to_h3k27ac_ratio/tCD34-high-top110/00_ologram_stats.tsv
        out/gtftk/merge_ologram_stats/ologram-h3k27ac-ratio-top110-classes.pdf

bed-cd34-ec-common-atac-peaks-classes-according-to-h3k27ac-top110
out/r/order_cd34_ec_common_atac_peaks_according_to_h3k27ac_ratio/{sample}", sample=['tCD34-high-top110.bed
    """
    input:
        stats = lambda wildcards: eval(mwconf['ids'][wildcards.ologram_id])
    output:
        pdf   = "out/{tool}{extra}/{ologram_id}.pdf"
    params:
        extra = params_extra
    threads:
        16
    wildcard_constraints:
        tool="gtftk/merge_ologram_stats"
    conda:
        "../envs/pygtftk.yaml"
    shell:
        "gtftk merge_ologram_stats {params.extra} -o {output.pdf} {input.stats}"

