rule deepTools_computeMatrix_extra:
    """
    Created:
        2018-10-29 12:03:28
    Aim:
        Can work both on reference-point and scale-regions
    Doc:
        https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html#reference-point
    Test:
        out/deepTools/plotHeatmap_sortRegions-descend_sortUsing-region_length_averageTypeSummaryPlot-mean_colorList-blueCyanYellowOrangeRed_heatmapHeight-28_heatmapWidth-1_whatToShow-h_xAxisLabel-peak-center_refPointLabel-0_boxAroundHeatmaps-no/
       out/deepTools/computeMatrix_reference-point_--referencePoint_center_-b_2000_-a_2000_-bs_200_--sortRegions_keep_-R_hg38-dClust-rowFeature-ealpha-and-20-random-regions_-S_hg38-H3K27ac-thymus-merged-wiq.txt.gz
       out/deepTools/plotHeatmap/deepTools/computeMatrix_reference-point_--referencePoint_center_-b_5000_-a_5000_--sortRegions_keep_-R_hg38-atac-peaks-from-wilson_-S_hg38-h3k27ac-cd34-ec-quantile-normalized.pdf
       out/deepTools/plotHeatmap/deepTools/compute_Matrix_scale-regions_bed-cd34-ec-common-atac-peaks-sorted-according-to-h3k27ac_bw-hg38-h3k27ac-cd34-ec-quantile-normalized.pdf
    """
    input:
        bw  = lambda wildcards: eval(config['ids'][wildcards.bw_list_id]),
        bed = lambda wildcards: eval(config['ids'][wildcards.bed_list_id])
    output:
        matrix = "out/{tool}{extra}_{bed_list_id}_{bw_list_id}.txt.gz"
    log:
                 "out/{tool}{extra}_{bed_list_id}_{bw_list_id}.log"
    benchmark:
                 "out/{tool}{extra}_{bed_list_id}_{bw_list_id}.benchmark.tsv"
    params:
        extra = params_extra
    conda:
        "../envs/deeptools.yaml"
    wildcard_constraints:
        tool="deepTools/computeMatrix",
        bw_list_id = "bw-[a-zA-Z0-9-]+",
        bed_list_id = "bed-[a-zA-Z0-9-]+"
    threads:
        16
    shell:
        "computeMatrix {params.extra} --regionsFileName {input.bed} "
        "--scoreFileName {input.bw} --outFileName {output.matrix} "
        "--numberOfProcessors {threads} &> {log}"

