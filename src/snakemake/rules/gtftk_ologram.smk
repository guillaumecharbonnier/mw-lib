rule gtftk_ologram:
    """
    Created:
        2019-04-02 16:19:18
    Aim:
        OverLap Of Genomic Regions Analysis using Monte Carlo to annotate peaks.
    Doc:
        https://dputhier.github.io/pygtftk/annotation.html#ologram
    Test:
        out/gtftk/ologram_tfbs_chrominfo-hg19_gtf-hg19-tfbs/bedtools/merge_-d_5/sort/_-k1,1_-k2,2n/crossmap/chain-hg38-to-hg19/cp/alias/experiments/hg38_atac_peaks_from_wilson/CD34/00_ologram_diagrams.pdf
        out/gtftk/ologram_tfbs_chrominfo-hg19_gtf-hg19-tfbs/bedtools/merge_-d_5/sort/_-k1,1_-k2,2n/crossmap/chain-hg38-to-hg19/r/dynamic_enhancers_in_thymopoiesis/dClust/rowFeature-all_cpg_hypo_meth_call__no_rmsk_mxy__no_donor_effect__distal/sortingMethod-kmeans_centers-8_nstart-100_itermax-200000_algorithm-Lloyd/sortingFeature-H3K27ac_peaks/subSortingMethod-kmeans_centers-32_nstart-100_itermax-200000_algorithm-Lloyd/subSortingFeature-ATAC_peaks/cluster-1/00_ologram_diagrams.pdf
        out/gtftk/ologram_tfbs_chrominfo-hg19_gtf-hg19-tfbs/bedtools/merge_-d_5/sort/_-k1,1_-k2,2n/crossmap/chain-hg38-to-hg19/r/dynamic_enhancers_in_thymopoiesis/dClust/rowFeature-no_rmsk_mxy__no_donor_effect__distal/sortingMethod-kmeans_centers-8_nstart-100_itermax-200000_algorithm-Lloyd/sortingFeature-cpg_meth_call/subSortingMethod-kmeans_centers-32_nstart-100_itermax-200000_algorithm-Lloyd/subSortingFeature-H3K27ac_peaks.ATAC_peaks/cluster-1/00_ologram_diagrams.pdf
    """
    input:
        chromInfo = lambda wildcards: config['ids'][wildcards.chrominfo_id],
        gtf = lambda wildcards: config['ids'][wildcards.gtf_id],
        bed="out/{filler}.bed"
    output:
        pdf       = "out/{tool}{extra}_{chrominfo_id}_{gtf_id}/{filler}/00_ologram_diagrams.pdf",
        stats     = "out/{tool}{extra}_{chrominfo_id}_{gtf_id}/{filler}/00_ologram_stats.tsv"
    params:
        extra = params_extra,
        outputdir = "out/{tool}{extra}_{chrominfo_id}_{gtf_id}/{filler}"
    threads:
        16
    wildcard_constraints:
        tool="gtftk/ologram",
        gtf_id="gtf-[a-zA-Z0-9-]+",
        chrominfo_id="chrominfo-[a-zA-Z0-9-]+"
    conda:
        "../envs/pygtftk.yaml"
    shell:
        """
        gtftk ologram\
            --inputfile {input.gtf} \
            --outputdir {params.outputdir} \
            --chrom-info {input.chromInfo} \
            --peak-file {input.bed} \
            --nb-threads {threads} \
            --no-date {params.extra}
        """
