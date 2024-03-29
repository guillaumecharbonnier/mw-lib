rule deepTools_multiBigwigSummary_BED_extra:
    """
    Created:
        2016-08-24 15h07
    Modified:
        2020-03-23 12:45:58 - adding eval() for input.bed. Will break all uses before this date.
    Doc:
        https://deeptools.readthedocs.io/en/develop/content/tools/multiBigwigSummary.html
    Note:
        input.bed can also be a gtf according to multiBigwigSummary documentation.
    Test:
        out/deepTools/multiBigwigSummary_BED-hg38-methylation-filtered-sites-in-thymus/Blueprint-thymic-populations-BS.npz
        out/deepTools/multiBigwigSummary_BED-hg38-atac-peaks-from-wilson-cd34-ec/hg38-h3k27ac-cd34-ec-quantile-normalized.tsv
        out/deepTools/multiBigwigSummary_BED_bed-hg38-ensembl-r95-ealpha/bw-hg38-tall-h3k27ac-qnorm.tsv
    """
    input:
        bw_list = lambda wildcards: eval(mwconf['ids'][wildcards.bw_list_id]),
        bed     = lambda wildcards: eval(mwconf['ids'][wildcards.bed_or_gtf_id])
    output:
        npz="out/{tool}{extra}_{bed_or_gtf_id}/{bw_list_id}.npz",
        tsv="out/{tool}{extra}_{bed_or_gtf_id}/{bw_list_id}.tsv"
    log:
            "out/{tool}{extra}_{bed_or_gtf_id}/{bw_list_id}.log"
    benchmark:
            "out/{tool}{extra}_{bed_or_gtf_id}/{bw_list_id}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="deepTools/multiBigwigSummary_BED"
    threads:
        MAX_THREADS
    conda:
        "../envs/deeptools.yaml"
    shell:
        "multiBigwigSummary BED-file --BED {input.bed} "
        "{params.extra} "
        "--bwfiles {input.bw_list} --outFileName {output.npz} "
        "--outRawCounts {output.tsv} --numberOfProcessors {threads} &> {log}"

rule deepTools_multiBigwigSummary_bins_extra:
    """
    Created:
        2016-08-24 15h07
    Aim:
        See if small structures are more correlated with nucleosomes from round spermatid or from condensing spermatid.
    Test:
        out/deepTools/multiBigwigSummary_bins/hg38-h3k27ac-cd34-ec-quantile-normalized.tsv
    """
    input:
        bw_list = lambda wildcards: eval(mwconf['ids'][wildcards.bw_list_id])
    output:
        npz="out/{tool}{extra}/{bw_list_id}.npz",
        tsv="out/{tool}{extra}/{bw_list_id}.tsv"
    log:
            "out/{tool}{extra}/{bw_list_id}.log"
    benchmark:
            "out/{tool}{extra}/{bw_list_id}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="deepTools/multiBigwigSummary_bins"
    threads:
        MAX_THREADS
    conda:
        "../envs/deeptools.yaml"
    shell:
        "multiBigwigSummary bins --bwfiles {input.bw_list} "
        "--outFileName {output.npz} --outRawCounts {output.tsv} "
        "--numberOfProcessors {threads} &> {log}"

