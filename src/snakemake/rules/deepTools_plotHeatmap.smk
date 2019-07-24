rule deepTools_plotHeatmap_extra:
    """
    Created:
        2018-10-29 13:37:50
    Test:
        out/deepTools/plotHeatmap_--sortRegions_keep_-colorList_blueCyanYellowOrangeRed_-whatToShow_phc_-xAxisLabel_peak-center_-refPointLabel_0/deepTools/computeMatrix_reference-point_--referencePoint_center_-b_2000_-a_2000_-bs_200_--sortRegions-keep_-R_hg38-dClust-rowFeature-ealpha-and-20-random-regions_-S_hg38-H3K27ac-thymus-merged-wiq.pdf
        out/deepTools/plotHeatmap_--sortRegions_keep_--samplesLabel_test_-colorList_blueCyanYellowOrangeRed_--whatToShow_phc_--xAxisLabel_peak-center_--refPointLabel_0/deepTools/computeMatrix_reference-point_--referencePoint_center_-b_2000_-a_2000_-bs_200_--sortRegions_keep_-R_hg38-dClust-rowFeature-ealpha-and-20-random-regions_-S_hg38-H3K27ac-thymus-merged-wiq.pdf
        out/deepTools/plotHeatmap_47aa/deepTools/computeMatrix_reference-point_--referencePoint_center_-b_2000_-a_2000_-bs_200_--sortRegions_keep_-R_hg38-dClust-rowFeature-ealpha-and-20-random-regions_-S_hg38-H3K27ac-thymus-merged-wiq.pdf
    Note:
        Force interpolationMethod as nearest for BS.
    """
    input:
        matrix="out/{filler}.txt.gz",
    output:
        pdf="out/{tool}{extra}/{filler}.pdf",
        bed="out/{tool}{extra}/{filler}.bed"
    log:
            "out/{tool}{extra}/{filler}.log"
    benchmark:
            "out/{tool}{extra}/{filler}.benchmark.tsv"
    threads:
        1
    wildcard_constraints:
        tool="deepTools/plotHeatmap"
    conda:
        "../envs/deeptools.yaml"
    priority:
        2
    params:
        extra = params_extra
    shell:
        "plotHeatmap --matrixFile {input.matrix} --outFileName {output.pdf} --outFileSortedRegions {output.bed} {params} &> {log}"

#####################
# ONLY LEGACY BELOW #
#####################

rule deepTools_plotHeatmap_sortRegions_sortUsing_averageTypeSummaryPlot_colorList_heatmapHeight_heatmapWidth_whatToShow_xAxisLabel_refPointLabel:
    """
    Created:
        2017-08-02 10:23:31
    Aim:
        Trying to find the settings that allow to plot dotted line for feature of interest.
    Test:
        out/deepTools/plotHeatmap_sortRegions-descend_sortUsing-region_length_averageTypeSummaryPlot-mean_colorList-blueCyanYellowOrangeRed_heatmapHeight-28_heatmapWidth-4_whatToShow-phc_xAxisLabel-tss_refPointLabel-0/deepTools/computeMatrix_reference-point_referencePoint-TSS_beforeRegionStartLength-1000_afterRegionStartLength-1000_bed-mm10-updownreg-ko-nut-r_bw-mm10-H4K5ac-Nut.pdf

        out/deepTools/plotHeatmap_sortRegions-descend_sortUsing-region_length_averageTypeSummaryPlot-mean_colorList-blueCyanYellowOrangeRed_heatmapHeight-28_heatmapWidth-4_whatToShow-phc_xAxisLabel-tss_refPointLabel-0/deepTools/computeMatrix_reference-point_referencePoint-TSS_beforeRegionStartLength-1000_afterRegionStartLength-1000_bed-mm10-classes-of-genes-by-rna-seq-expression_bw-mm10-MNase-wt.pdf

        out/deepTools/plotHeatmap_sortRegions-descend_sortUsing-region_length_averageTypeSummaryPlot-mean_colorList-blueCyanYellowOrangeRed_heatmapHeight-28_heatmapWidth-4_whatToShow-phc_xAxisLabel-tss_refPointLabel-0/deepTools/computeMatrix_reference-point_referencePoint-center_beforeRegionStartLength-2000_afterRegionStartLength-2000_bed-hg38-macs2-peaks-H3K27ac-thymus_bw-hg38-H3K27ac-thymus.pdf
    """
    input:
        matrix="out/{filler}.txt.gz",
        plotHeatmap="opt/miniconda/envs/deeptools/bin/plotHeatmap"
    output:
        pdf="out/deepTools/plotHeatmap_sortRegions-{deepTools_plotHeatmap_sortRegions}_sortUsing-{deepTools_plotHeatmap_sortUsing}_averageTypeSummaryPlot-{deepTools_plotHeatmap_averageTypeSummaryPlot}_colorList-{deepTools_plotHeatmap_colorList_id}_heatmapHeight-{deepTools_plotHeatmap_heatmapHeight}_heatmapWidth-{deepTools_plotHeatmap_heatmapWidth}_whatToShow-{deepTools_plotHeatmap_whatToShow_id}_xAxisLabel-{deepTools_plotHeatmap_xAxisLabel_id}_refPointLabel-{deepTools_plotHeatmap_refPointLabel_id}/{filler}.pdf",
        bed="out/deepTools/plotHeatmap_sortRegions-{deepTools_plotHeatmap_sortRegions}_sortUsing-{deepTools_plotHeatmap_sortUsing}_averageTypeSummaryPlot-{deepTools_plotHeatmap_averageTypeSummaryPlot}_colorList-{deepTools_plotHeatmap_colorList_id}_heatmapHeight-{deepTools_plotHeatmap_heatmapHeight}_heatmapWidth-{deepTools_plotHeatmap_heatmapWidth}_whatToShow-{deepTools_plotHeatmap_whatToShow_id}_xAxisLabel-{deepTools_plotHeatmap_xAxisLabel_id}_refPointLabel-{deepTools_plotHeatmap_refPointLabel_id}/{filler}.bed"
    threads:
        1
    priority:
        2
    params:
        deepTools_plotHeatmap_colorList              = params_deepTools_plotHeatmap_colorList,
        deepTools_plotHeatmap_whatToShow             = params_deepTools_plotHeatmap_whatToShow,
        deepTools_plotHeatmap_xAxisLabel             = params_deepTools_plotHeatmap_xAxisLabel,
        deepTools_plotHeatmap_refPointLabel          = params_deepTools_plotHeatmap_refPointLabel,
        deepTools_plotHeatmap_sortRegions            = "--sortRegion {deepTools_plotHeatmap_sortRegions}",
        deepTools_plotHeatmap_sortUsing              = "--sortUsing {deepTools_plotHeatmap_sortUsing}",
        deepTools_plotHeatmap_averageTypeSummaryPlot = "--averageTypeSummaryPlot {deepTools_plotHeatmap_averageTypeSummaryPlot}",
        deepTools_plotHeatmap_heatmapHeight          = "--heatmapHeight {deepTools_plotHeatmap_heatmapHeight}",
        deepTools_plotHeatmap_heatmapWidth           = "--heatmapWidth {deepTools_plotHeatmap_heatmapWidth}",
    shell:
        "{input.plotHeatmap} --matrixFile {input.matrix} --outFileName {output.pdf} --outFileSortedRegions {output.bed} {params}"

rule deepTools_plotHeatmap_sortRegions_sortUsing_averageTypeSummaryPlot_missingDataColor_colorList_heatmapHeight_heatmapWidth_whatToShow_xAxisLabel_refPointLabel_boxAroundHeatmaps:
    """
    Created:
        2018-04-09 15:42:46
    Aim:
        Adding the missingDataColor parameter to see if Bisulfite-Seq data is cleaner with white as missigDataColor instead of black.
    Test:
    """
    input:
        matrix="out/{filler}.txt.gz",
        plotHeatmap="opt/miniconda/envs/deeptools/bin/plotHeatmap"
    output:
        pdf="out/deepTools/plotHeatmap_sortRegions-{deepTools_plotHeatmap_sortRegions}_sortUsing-{deepTools_plotHeatmap_sortUsing}_averageTypeSummaryPlot-{deepTools_plotHeatmap_averageTypeSummaryPlot}_missingDataColor-{deepTools_plotHeatmap_missingDataColor}_colorList-{deepTools_plotHeatmap_colorList_id}_heatmapHeight-{deepTools_plotHeatmap_heatmapHeight}_heatmapWidth-{deepTools_plotHeatmap_heatmapWidth}_whatToShow-{deepTools_plotHeatmap_whatToShow_id}_xAxisLabel-{deepTools_plotHeatmap_xAxisLabel_id}_refPointLabel-{deepTools_plotHeatmap_refPointLabel_id}_boxAroundHeatmaps-{deepTools_plotHeatmap_boxAroundHeatmaps}/{filler}.pdf",
        bed="out/deepTools/plotHeatmap_sortRegions-{deepTools_plotHeatmap_sortRegions}_sortUsing-{deepTools_plotHeatmap_sortUsing}_averageTypeSummaryPlot-{deepTools_plotHeatmap_averageTypeSummaryPlot}_missingDataColor-{deepTools_plotHeatmap_missingDataColor}_colorList-{deepTools_plotHeatmap_colorList_id}_heatmapHeight-{deepTools_plotHeatmap_heatmapHeight}_heatmapWidth-{deepTools_plotHeatmap_heatmapWidth}_whatToShow-{deepTools_plotHeatmap_whatToShow_id}_xAxisLabel-{deepTools_plotHeatmap_xAxisLabel_id}_refPointLabel-{deepTools_plotHeatmap_refPointLabel_id}_boxAroundHeatmaps-{deepTools_plotHeatmap_boxAroundHeatmaps}/{filler}.bed"
    threads:
        1
    priority:
        2
    params:
        deepTools_plotHeatmap_colorList              = params_deepTools_plotHeatmap_colorList,
        deepTools_plotHeatmap_whatToShow             = params_deepTools_plotHeatmap_whatToShow,
        deepTools_plotHeatmap_xAxisLabel             = params_deepTools_plotHeatmap_xAxisLabel,
        deepTools_plotHeatmap_refPointLabel          = params_deepTools_plotHeatmap_refPointLabel,
        deepTools_plotHeatmap_sortRegions            = "--sortRegion {deepTools_plotHeatmap_sortRegions}",
        deepTools_plotHeatmap_sortUsing              = "--sortUsing {deepTools_plotHeatmap_sortUsing}",
        deepTools_plotHeatmap_averageTypeSummaryPlot = "--averageTypeSummaryPlot {deepTools_plotHeatmap_averageTypeSummaryPlot}",
        deepTools_plotHeatmap_missingDataColor       = "--missingDataColor {deepTools_plotHeatmap_missingDataColor}",
        deepTools_plotHeatmap_heatmapHeight          = "--heatmapHeight {deepTools_plotHeatmap_heatmapHeight}",
        deepTools_plotHeatmap_heatmapWidth           = "--heatmapWidth {deepTools_plotHeatmap_heatmapWidth}",
        deepTools_plotHeatmap_boxAroundHeatmaps      = "--boxAroundHeatmaps {deepTools_plotHeatmap_boxAroundHeatmaps}"
    shell:
        "{input.plotHeatmap} --matrixFile {input.matrix} --outFileName {output.pdf} --outFileSortedRegions {output.bed} {params}"

rule deepTools_plotHeatmap_sortRegions_sortUsing_averageTypeSummaryPlot_colorList_heatmapHeight_heatmapWidth_whatToShow_xAxisLabel_refPointLabel_boxAroundHeatmaps:
    """
    Created:
        2017-08-02 10:23:31
    Aim:
        Trying to find the settings that allow to plot dotted line for feature of interest.
    Test:
    """
    input:
        matrix="out/{filler}.txt.gz",
        plotHeatmap="opt/miniconda/envs/deeptools/bin/plotHeatmap"
    output:
        pdf="out/deepTools/plotHeatmap_sortRegions-{deepTools_plotHeatmap_sortRegions}_sortUsing-{deepTools_plotHeatmap_sortUsing}_averageTypeSummaryPlot-{deepTools_plotHeatmap_averageTypeSummaryPlot}_colorList-{deepTools_plotHeatmap_colorList_id}_heatmapHeight-{deepTools_plotHeatmap_heatmapHeight}_heatmapWidth-{deepTools_plotHeatmap_heatmapWidth}_whatToShow-{deepTools_plotHeatmap_whatToShow_id}_xAxisLabel-{deepTools_plotHeatmap_xAxisLabel_id}_refPointLabel-{deepTools_plotHeatmap_refPointLabel_id}_boxAroundHeatmaps-{deepTools_plotHeatmap_boxAroundHeatmaps}/{filler}.pdf",
        bed="out/deepTools/plotHeatmap_sortRegions-{deepTools_plotHeatmap_sortRegions}_sortUsing-{deepTools_plotHeatmap_sortUsing}_averageTypeSummaryPlot-{deepTools_plotHeatmap_averageTypeSummaryPlot}_colorList-{deepTools_plotHeatmap_colorList_id}_heatmapHeight-{deepTools_plotHeatmap_heatmapHeight}_heatmapWidth-{deepTools_plotHeatmap_heatmapWidth}_whatToShow-{deepTools_plotHeatmap_whatToShow_id}_xAxisLabel-{deepTools_plotHeatmap_xAxisLabel_id}_refPointLabel-{deepTools_plotHeatmap_refPointLabel_id}_boxAroundHeatmaps-{deepTools_plotHeatmap_boxAroundHeatmaps}/{filler}.bed"
    threads:
        1
    priority:
        2
    params:
        deepTools_plotHeatmap_colorList              = params_deepTools_plotHeatmap_colorList,
        deepTools_plotHeatmap_whatToShow             = params_deepTools_plotHeatmap_whatToShow,
        deepTools_plotHeatmap_xAxisLabel             = params_deepTools_plotHeatmap_xAxisLabel,
        deepTools_plotHeatmap_refPointLabel          = params_deepTools_plotHeatmap_refPointLabel,
        deepTools_plotHeatmap_sortRegions            = "--sortRegion {deepTools_plotHeatmap_sortRegions}",
        deepTools_plotHeatmap_sortUsing              = "--sortUsing {deepTools_plotHeatmap_sortUsing}",
        deepTools_plotHeatmap_averageTypeSummaryPlot = "--averageTypeSummaryPlot {deepTools_plotHeatmap_averageTypeSummaryPlot}",
        deepTools_plotHeatmap_heatmapHeight          = "--heatmapHeight {deepTools_plotHeatmap_heatmapHeight}",
        deepTools_plotHeatmap_heatmapWidth           = "--heatmapWidth {deepTools_plotHeatmap_heatmapWidth}",
        deepTools_plotHeatmap_boxAroundHeatmaps      = "--boxAroundHeatmaps {deepTools_plotHeatmap_boxAroundHeatmaps}"
    shell:
        "{input.plotHeatmap} --matrixFile {input.matrix} --outFileName {output.pdf} --outFileSortedRegions {output.bed} {params}"

rule deepTools_plotHeatmap_kmeans_sortRegions_sortUsing_averageTypeSummaryPlot_colorList_heatmapHeight_heatmapWidth_whatToShow_xAxisLabel_refPointLabel:
    """
    Created:
        2017-05-08 16:37:35
    Test:
        out/deepTools/plotHeatmap_kmeans-0_sortRegions-descend_sortUsing-mean_averageTypeSummaryPlot-mean_colorList-blueCyanYellowOrangeRed_heatmapHeight-28_heatmapWidth-4_whatToShow-phc_xAxisLabel-peak-center_refPointLabel-0/deepTools/computeMatrix_reference-point_referencePoint-center_beforeRegionStartLength-1000_afterRegionStartLength-1000_bed-mm10-h2alap1-peaks_bw-mm10-MNase-wt.pdf
    """
    input:
        matrix="out/{filler}.txt.gz",
        plotHeatmap="opt/miniconda/envs/deeptools/bin/plotHeatmap"
    output:
        pdf="out/deepTools/plotHeatmap_kmeans-{deepTools_plotHeatmap_kmeans}_sortRegions-{deepTools_plotHeatmap_sortRegions}_sortUsing-{deepTools_plotHeatmap_sortUsing}_averageTypeSummaryPlot-{deepTools_plotHeatmap_averageTypeSummaryPlot}_colorList-{deepTools_plotHeatmap_colorList_id}_heatmapHeight-{deepTools_plotHeatmap_heatmapHeight}_heatmapWidth-{deepTools_plotHeatmap_heatmapWidth}_whatToShow-{deepTools_plotHeatmap_whatToShow_id}_xAxisLabel-{deepTools_plotHeatmap_xAxisLabel_id}_refPointLabel-{deepTools_plotHeatmap_refPointLabel_id}/{filler}.pdf",
        bed="out/deepTools/plotHeatmap_kmeans-{deepTools_plotHeatmap_kmeans}_sortRegions-{deepTools_plotHeatmap_sortRegions}_sortUsing-{deepTools_plotHeatmap_sortUsing}_averageTypeSummaryPlot-{deepTools_plotHeatmap_averageTypeSummaryPlot}_colorList-{deepTools_plotHeatmap_colorList_id}_heatmapHeight-{deepTools_plotHeatmap_heatmapHeight}_heatmapWidth-{deepTools_plotHeatmap_heatmapWidth}_whatToShow-{deepTools_plotHeatmap_whatToShow_id}_xAxisLabel-{deepTools_plotHeatmap_xAxisLabel_id}_refPointLabel-{deepTools_plotHeatmap_refPointLabel_id}/{filler}.bed"
    threads:
        1
    priority:
        2
    params:
        deepTools_plotHeatmap_colorList              = params_deepTools_plotHeatmap_colorList,
        deepTools_plotHeatmap_whatToShow             = params_deepTools_plotHeatmap_whatToShow,
        deepTools_plotHeatmap_xAxisLabel             = params_deepTools_plotHeatmap_xAxisLabel,
        deepTools_plotHeatmap_refPointLabel          = params_deepTools_plotHeatmap_refPointLabel,
        deepTools_plotHeatmap_kmeans                 = "--kmeans {deepTools_plotHeatmap_kmeans}",
        deepTools_plotHeatmap_sortRegions            = "--sortRegion {deepTools_plotHeatmap_sortRegions}",
        deepTools_plotHeatmap_sortUsing              = "--sortUsing {deepTools_plotHeatmap_sortUsing}",
        deepTools_plotHeatmap_averageTypeSummaryPlot = "--averageTypeSummaryPlot {deepTools_plotHeatmap_averageTypeSummaryPlot}",
        deepTools_plotHeatmap_heatmapHeight          = "--heatmapHeight {deepTools_plotHeatmap_heatmapHeight}",
        deepTools_plotHeatmap_heatmapWidth           = "--heatmapWidth {deepTools_plotHeatmap_heatmapWidth}",
    shell:
        "{input.plotHeatmap} --matrixFile {input.matrix} --outFileName {output.pdf} --outFileSortedRegions {output.bed} {params}"

rule deepTools_plotHeatmap_sortRegions_sortUsing_averageTypeSummaryPlot_colorList_zMin_zMax_heatmapHeight_heatmapWidth_whatToShow_xAxisLabel_refPointLabel_boxAroundHeatmaps:
    """
    Created:
        2018-03-21 18:46:27
    Aim:
        Trying to find the settings that allow to plot dotted line for feature of interest.
    Test:
    """
    input:
        matrix="out/{filler}.txt.gz",
        plotHeatmap="opt/miniconda/envs/deeptools/bin/plotHeatmap"
    output:
        pdf="out/deepTools/plotHeatmap_sortRegions-{deepTools_plotHeatmap_sortRegions}_sortUsing-{deepTools_plotHeatmap_sortUsing}_averageTypeSummaryPlot-{deepTools_plotHeatmap_averageTypeSummaryPlot}_colorList-{deepTools_plotHeatmap_colorList_id}_zMin-{deepTools_plotHeatmap_zMin}_zMax-{deepTools_plotHeatmap_zMax}_heatmapHeight-{deepTools_plotHeatmap_heatmapHeight}_heatmapWidth-{deepTools_plotHeatmap_heatmapWidth}_whatToShow-{deepTools_plotHeatmap_whatToShow_id}_xAxisLabel-{deepTools_plotHeatmap_xAxisLabel_id}_refPointLabel-{deepTools_plotHeatmap_refPointLabel_id}_boxAroundHeatmaps-{deepTools_plotHeatmap_boxAroundHeatmaps}/{filler}.pdf",
        bed="out/deepTools/plotHeatmap_sortRegions-{deepTools_plotHeatmap_sortRegions}_sortUsing-{deepTools_plotHeatmap_sortUsing}_averageTypeSummaryPlot-{deepTools_plotHeatmap_averageTypeSummaryPlot}_colorList-{deepTools_plotHeatmap_colorList_id}_zMin-{deepTools_plotHeatmap_zMin}_zMax-{deepTools_plotHeatmap_zMax}_heatmapHeight-{deepTools_plotHeatmap_heatmapHeight}_heatmapWidth-{deepTools_plotHeatmap_heatmapWidth}_whatToShow-{deepTools_plotHeatmap_whatToShow_id}_xAxisLabel-{deepTools_plotHeatmap_xAxisLabel_id}_refPointLabel-{deepTools_plotHeatmap_refPointLabel_id}_boxAroundHeatmaps-{deepTools_plotHeatmap_boxAroundHeatmaps}/{filler}.bed"
    threads:
        1
    priority:
        2
    params:
        deepTools_plotHeatmap_colorList              = params_deepTools_plotHeatmap_colorList,
        deepTools_plotHeatmap_whatToShow             = params_deepTools_plotHeatmap_whatToShow,
        deepTools_plotHeatmap_xAxisLabel             = params_deepTools_plotHeatmap_xAxisLabel,
        deepTools_plotHeatmap_refPointLabel          = params_deepTools_plotHeatmap_refPointLabel,
        deepTools_plotHeatmap_sortRegions            = "--sortRegion {deepTools_plotHeatmap_sortRegions}",
        deepTools_plotHeatmap_sortUsing              = "--sortUsing {deepTools_plotHeatmap_sortUsing}",
        deepTools_plotHeatmap_averageTypeSummaryPlot = "--averageTypeSummaryPlot {deepTools_plotHeatmap_averageTypeSummaryPlot}",
        deepTools_plotHeatmap_zMin                   = "--zMin {deepTools_plotHeatmap_zMin}",
        deepTools_plotHeatmap_zMax                   = "--zMax {deepTools_plotHeatmap_zMax}",
        deepTools_plotHeatmap_heatmapHeight          = "--heatmapHeight {deepTools_plotHeatmap_heatmapHeight}",
        deepTools_plotHeatmap_heatmapWidth           = "--heatmapWidth {deepTools_plotHeatmap_heatmapWidth}",
        deepTools_plotHeatmap_boxAroundHeatmaps      = "--boxAroundHeatmaps {deepTools_plotHeatmap_boxAroundHeatmaps}"
    shell:
        "{input.plotHeatmap} --matrixFile {input.matrix} --outFileName {output.pdf} --outFileSortedRegions {output.bed} {params}"


####
# IMPORTANT : rules below this point should not work anymore since 2018-04-09 17:29:40.
# To make them work, update rule skeleton according to rules above.
# I am too lazy to make the update myself if I do not need them right now.
####

rule deepTools_plotHeatmap_sortRegions_sortUsing_averageTypeSummaryPlot_colorList_heatmapHeight_heatmapWidth_whatToShow_xAxisLabel_refPointLabel_perGroup:
    """
    Created:
        2017-10-24 10:07:02
    Aim:
        PerGroup is nice in some cases.
    Test:
    """
    input:
        matrix="out/{filler}.txt.gz",
        plotHeatmap="opt/miniconda/envs/deeptools/bin/plotHeatmap"
    output:
        pdf="out/deepTools/plotHeatmap_sortRegions-{sortRegions}_sortUsing-{sortUsing}_averageTypeSummaryPlot-{averageTypeSummaryPlot}_colorList-{colorList_id}_heatmapHeight-{heatmapHeight}_heatmapWidth-{heatmapWidth}_whatToShow-{whatToShow_id}_xAxisLabel-{xAxisLabel_id}_refPointLabel-{refPointLabel_id}_perGroup/{filler}.pdf",
        bed="out/deepTools/plotHeatmap_sortRegions-{sortRegions}_sortUsing-{sortUsing}_averageTypeSummaryPlot-{averageTypeSummaryPlot}_colorList-{colorList_id}_heatmapHeight-{heatmapHeight}_heatmapWidth-{heatmapWidth}_whatToShow-{whatToShow_id}_xAxisLabel-{xAxisLabel_id}_refPointLabel-{refPointLabel_id}_perGroup/{filler}.bed"
    threads:
        1
    priority:
        2
    params:
        colorList=params_deepTools_plotHeatmap_colorList,
        whatToShow=params_deepTools_plotHeatmap_whatToShow,
        xAxisLabel=params_deepTools_plotHeatmap_xAxisLabel,
        refPointLabel=params_deepTools_plotHeatmap_refPointLabel,
    wildcard_constraints:
        sortRegions="descend|ascend|no",
        sortUsing="mean|median|max|min|sum|region_length",
        averageTypeSummaryPlot="mean|median|min|max|sum|std",
        colorList_id="[a-zA-Z]+",
        heatmapHeight="[0-9]+", #Plot height in cm. The default value is 28. The minimum value is 3 and the maximum is 100.
        heatmapWidth="[0-9]+", #Plot width in cm. The default value is 4 The minimum value is 1 and the maximum is 100.
        whatToShow_id="phc|ph|h|hc", # write a function to map these to path compliant ids.
        xAxisLabel_id="[a-zA-Z-]+", # write a function to map these to path compliant ids.
        refPointLabel_id="[a-zA-Z0-9-]+", #  write a function to map these to path compliant ids.
    shell:
        """
        {input.plotHeatmap} \
            --matrixFile {input.matrix} \
            --outFileName {output.pdf} \
            --outFileSortedRegions {output.bed} \
            --sortRegions {wildcards.sortRegions} \
            --sortUsing {wildcards.sortUsing} \
            --averageTypeSummaryPlot {wildcards.averageTypeSummaryPlot} \
            --colorList {params.colorList} \
            --heatmapHeight {wildcards.heatmapHeight} \
            --heatmapWidth {wildcards.heatmapWidth} \
            --whatToShow {params.whatToShow} \
            --xAxisLabel {params.xAxisLabel} \
            --refPointLabel {params.refPointLabel} \
            --perGroup
        """

rule deepTools_plotHeatmap_kmeans_sortRegions_sortUsing_averageTypeSummaryPlot_colorList_heatmapHeight_heatmapWidth_whatToShow_xAxisLabel_refPointLabel_perGroup:
    """
    Created:
        2017-10-24 10:18:35
    Aim:
        PerGroup is nice in some cases.
    Test:
    """
    input:
        matrix="out/{filler}.txt.gz",
        plotHeatmap="opt/miniconda/envs/deeptools/bin/plotHeatmap"
    output:
        pdf="out/deepTools/plotHeatmap_kmeans-{kmeans}_sortRegions-{sortRegions}_sortUsing-{sortUsing}_averageTypeSummaryPlot-{averageTypeSummaryPlot}_colorList-{colorList_id}_heatmapHeight-{heatmapHeight}_heatmapWidth-{heatmapWidth}_whatToShow-{whatToShow_id}_xAxisLabel-{xAxisLabel_id}_refPointLabel-{refPointLabel_id}_perGroup/{filler}.pdf",
        bed="out/deepTools/plotHeatmap_kmeans-{kmeans}_sortRegions-{sortRegions}_sortUsing-{sortUsing}_averageTypeSummaryPlot-{averageTypeSummaryPlot}_colorList-{colorList_id}_heatmapHeight-{heatmapHeight}_heatmapWidth-{heatmapWidth}_whatToShow-{whatToShow_id}_xAxisLabel-{xAxisLabel_id}_refPointLabel-{refPointLabel_id}_perGroup/{filler}.bed"
    threads:
        1
    priority:
        2
    params:
        colorList=params_deepTools_plotHeatmap_colorList,
        whatToShow=params_deepTools_plotHeatmap_whatToShow,
        xAxisLabel=params_deepTools_plotHeatmap_xAxisLabel,
        refPointLabel=params_deepTools_plotHeatmap_refPointLabel,
    wildcard_constraints:
        kmeans="[0-9]+",
        sortRegions="descend|ascend|no",
        sortUsing="mean|median|max|min|sum|region_length",
        averageTypeSummaryPlot="mean|median|min|max|sum|std",
        colorList_id="[a-zA-Z]+",
        heatmapHeight="[0-9]+", #Plot height in cm. The default value is 28. The minimum value is 3 and the maximum is 100.
        heatmapWidth="[0-9]+", #Plot width in cm. The default value is 4 The minimum value is 1 and the maximum is 100.
        whatToShow_id="phc|ph|h|hc", # write a function to map these to path compliant ids.
        xAxisLabel_id="[a-zA-Z-]+", # write a function to map these to path compliant ids.
        refPointLabel_id="[a-zA-Z0-9-]+", #  write a function to map these to path compliant ids.
    shell:
        """
        {input.plotHeatmap} \
            --matrixFile {input.matrix} \
            --outFileName {output.pdf} \
            --outFileSortedRegions {output.bed} \
            --kmeans {wildcards.kmeans} \
            --sortRegions {wildcards.sortRegions} \
            --sortUsing {wildcards.sortUsing} \
            --averageTypeSummaryPlot {wildcards.averageTypeSummaryPlot} \
            --colorList {params.colorList} \
            --heatmapHeight {wildcards.heatmapHeight} \
            --heatmapWidth {wildcards.heatmapWidth} \
            --whatToShow {params.whatToShow} \
            --xAxisLabel {params.xAxisLabel} \
            --refPointLabel {params.refPointLabel} \
            --perGroup
        """

rule deepTools_plotHeatmap_sortRegions_sortUsing_averageTypeSummaryPlot_colorList_heatmapHeight_heatmapWidth_whatToShow_xAxisLabel_refPointLabel_yMin_yMax:
    """
    Created:
        2017-10-09 11:33:01
    Aim:
        Checking strange scales on profile: the different samples seems to have different scales.
    Test:
        out/deepTools/plotHeatmap_sortRegions-descend_sortUsing-mean_averageTypeSummaryPlot-mean_colorList-blueCyanYellowOrangeRed_heatmapHeight-9_heatmapWidth-3_whatToShow-phc_xAxisLabel-tss_refPointLabel-0_yMin-0_yMax-100/deepTools/computeMatrix_reference-point_referencePoint-TSS_beforeRegionStartLength-3000_afterRegionStartLength-3000_bed-mm10-microarray-upreg-ko-nut_bw-mm10-H4K5ac-Nut.pdf
    """
    input:
        matrix="out/{filler}.txt.gz",
        plotHeatmap="opt/miniconda/envs/deeptools/bin/plotHeatmap"
    output:
        pdf="out/deepTools/plotHeatmap_sortRegions-{sortRegions}_sortUsing-{sortUsing}_averageTypeSummaryPlot-{averageTypeSummaryPlot}_colorList-{colorList_id}_heatmapHeight-{heatmapHeight}_heatmapWidth-{heatmapWidth}_whatToShow-{whatToShow_id}_xAxisLabel-{xAxisLabel_id}_refPointLabel-{refPointLabel_id}_yMin-{yMin}_yMax-{yMax}/{filler}.pdf",
        bed="out/deepTools/plotHeatmap_sortRegions-{sortRegions}_sortUsing-{sortUsing}_averageTypeSummaryPlot-{averageTypeSummaryPlot}_colorList-{colorList_id}_heatmapHeight-{heatmapHeight}_heatmapWidth-{heatmapWidth}_whatToShow-{whatToShow_id}_xAxisLabel-{xAxisLabel_id}_refPointLabel-{refPointLabel_id}_yMin-{yMin}_yMax-{yMax}/{filler}.bed"
    threads:
        1
    priority:
        2
    params:
        colorList=params_deepTools_plotHeatmap_colorList,
        whatToShow=params_deepTools_plotHeatmap_whatToShow,
        xAxisLabel=params_deepTools_plotHeatmap_xAxisLabel,
        refPointLabel=params_deepTools_plotHeatmap_refPointLabel,
    wildcard_constraints:
        sortRegions="descend|ascend|no",
        sortUsing="mean|median|max|min|sum|region_length",
        averageTypeSummaryPlot="mean|median|min|max|sum|std",
        colorList_id="[a-zA-Z]+",
        heatmapHeight="[0-9]+", #Plot height in cm. The default value is 28. The minimum value is 3 and the maximum is 100.
        heatmapWidth="[0-9]+", #Plot width in cm. The default value is 4 The minimum value is 1 and the maximum is 100.
        whatToShow_id="phc|ph|h|hc", # write a function to map these to path compliant ids.
        xAxisLabel_id="[a-zA-Z-]+", # write a function to map these to path compliant ids.
        refPointLabel_id="[a-zA-Z0-9-]+", #  write a function to map these to path compliant ids.
        yMin="[0-9-]+",
        yMax="[0-9-]+"
    shell:
        """
        {input.plotHeatmap} \
            --matrixFile {input.matrix} \
            --outFileName {output.pdf} \
            --outFileSortedRegions {output.bed} \
            --sortRegions {wildcards.sortRegions} \
            --sortUsing {wildcards.sortUsing} \
            --averageTypeSummaryPlot {wildcards.averageTypeSummaryPlot} \
            --colorList {params.colorList} \
            --heatmapHeight {wildcards.heatmapHeight} \
            --heatmapWidth {wildcards.heatmapWidth} \
            --whatToShow {params.whatToShow} \
            --xAxisLabel {params.xAxisLabel} \
            --refPointLabel {params.refPointLabel} \
            --yMin {wildcards.yMin} \
            --yMax {wildcards.yMax}
        """

rule deepTools_plotHeatmap_sortRegions_sortUsing_averageTypeSummaryPlot_colorList_zMin_zMax_heatmapHeight_heatmapWidth_whatToShow_xAxisLabel_refPointLabel_yMin_yMax:
    """
    Created:
        2017-10-10 12:54:52
    Aim:
        I need to center on 0 y and z scales
    Test:
        out/deepTools/plotHeatmap_sortRegions-descend_sortUsing-mean_averageTypeSummaryPlot-mean_colorList-blueCyanYellowOrangeRed_zMin--0.5_zMax-0.5_heatmapHeight-15_heatmapWidth-3_whatToShow-phc_xAxisLabel-tss_refPointLabel-0_yMin--0.5_yMax-0.5/deepTools/computeMatrix_reference-point_referencePoint-TSS_beforeRegionStartLength-3000_afterRegionStartLength-3000_bed-mm10-microarray-downreg-ko-nut/deepTools/bamCompare_scaleFactorsMethod-readCount_ratio-log2_binSize-10_extendReads-150/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run140_run141/H4K5ac-Nut-WT_vs_H4K5ac-Nut-KO.pdf
    Note:
        deepTools_plotHeatmap_sortRegions_sortUsing_averageTypeSummaryPlot_colorList_heatmapHeight_heatmapWidth_whatToShow_xAxisLabel_refPointLabel_yMin_yMax
        and
        deepTools_plotHeatmap_sortRegions_sortUsing_averageTypeSummaryPlot_colorList_zMin_zMax_heatmapHeight_heatmapWidth_whatToShow_xAxisLabel_refPointLabel_yMin_yMax are ambiguous for the file out/deepTools/plotHeatmap_sortRegions-descend_sortUsing-mean_averageTypeSummaryPlot-mean_colorList-blueCyanYellowOrangeRed_zMin--2_zMax-5_heatmapHeight-15_heatmapWidth-3_whatToShow-phc_xAxisLabel-tss_refPointLabel-0_yMin-1_yMax-3/deepTools/computeMatrix_reference-point_referencePoint-TSS_beforeRegionStartLength-3000_afterRegionStartLength-3000_bed-mm10-microarray-downreg-ko-nut/deepTools/bamCompare_scaleFactorsMethod-SES_ratio-log2_binSize-10_extendReads-150/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run140_run141/H4K5ac-Nut-KO_vs_Input-Nut-KO.pdf.

    """
    input:
        matrix="out/{filler}.txt.gz",
        plotHeatmap="opt/miniconda/envs/deeptools/bin/plotHeatmap"
    output:
        pdf="out/deepTools/plotHeatmap_sortRegions-{sortRegions}_sortUsing-{sortUsing}_averageTypeSummaryPlot-{averageTypeSummaryPlot}_colorList-{colorList_id}_zMin-{zMin}_zMax-{zMax}_heatmapHeight-{heatmapHeight}_heatmapWidth-{heatmapWidth}_whatToShow-{whatToShow_id}_xAxisLabel-{xAxisLabel_id}_refPointLabel-{refPointLabel_id}_yMin-{yMin}_yMax-{yMax}/{filler}.pdf",
        bed="out/deepTools/plotHeatmap_sortRegions-{sortRegions}_sortUsing-{sortUsing}_averageTypeSummaryPlot-{averageTypeSummaryPlot}_colorList-{colorList_id}_zMin-{zMin}_zMax-{zMax}_heatmapHeight-{heatmapHeight}_heatmapWidth-{heatmapWidth}_whatToShow-{whatToShow_id}_xAxisLabel-{xAxisLabel_id}_refPointLabel-{refPointLabel_id}_yMin-{yMin}_yMax-{yMax}/{filler}.bed"
    threads:
        1
    priority:
        2
    params:
        colorList=params_deepTools_plotHeatmap_colorList,
        whatToShow=params_deepTools_plotHeatmap_whatToShow,
        xAxisLabel=params_deepTools_plotHeatmap_xAxisLabel,
        refPointLabel=params_deepTools_plotHeatmap_refPointLabel,
    wildcard_constraints:
        sortRegions="descend|ascend|no",
        sortUsing="mean|median|max|min|sum|region_length",
        averageTypeSummaryPlot="mean|median|min|max|sum|std",
        colorList_id="[a-zA-Z]+",
        heatmapHeight="[0-9]+", #Plot height in cm. The default value is 28. The minimum value is 3 and the maximum is 100.
        heatmapWidth="[0-9]+", #Plot width in cm. The default value is 4 The minimum value is 1 and the maximum is 100.
        whatToShow_id="phc|ph|h|hc", # write a function to map these to path compliant ids.
        xAxisLabel_id="[a-zA-Z-]+", # write a function to map these to path compliant ids.
        refPointLabel_id="[a-zA-Z0-9-]+", #  write a function to map these to path compliant ids.
        zMin="[0-9-.]+",
        zMax="[0-9-.]+",
        yMin="[0-9-.]+",
        yMax="[0-9-.]+"
    shell:
        """
        {input.plotHeatmap} \
            --matrixFile {input.matrix} \
            --outFileName {output.pdf} \
            --outFileSortedRegions {output.bed} \
            --sortRegions {wildcards.sortRegions} \
            --sortUsing {wildcards.sortUsing} \
            --averageTypeSummaryPlot {wildcards.averageTypeSummaryPlot} \
            --colorList {params.colorList} \
            --heatmapHeight {wildcards.heatmapHeight} \
            --heatmapWidth {wildcards.heatmapWidth} \
            --whatToShow {params.whatToShow} \
            --xAxisLabel {params.xAxisLabel} \
            --refPointLabel {params.refPointLabel} \
            --zMin {wildcards.zMin} \
            --zMax {wildcards.zMax} \
            --yMin {wildcards.yMin} \
            --yMax {wildcards.yMax}
        """

rule deepTools_plotHeatmap_sortRegions_sortUsing_averageTypeSummaryPlot_colorList_heatmapHeight_heatmapWidth_whatToShow_xAxisLabel_startLabel_endLabel:
    """
    Created:
        2017-10-18 15:42:26
    Test:
    """
    input:
        matrix="out/{filler}.txt.gz",
        plotHeatmap="opt/miniconda/envs/deeptools/bin/plotHeatmap"
    output:
        pdf="out/deepTools/plotHeatmap_sortRegions-{sortRegions}_sortUsing-{sortUsing}_averageTypeSummaryPlot-{averageTypeSummaryPlot}_colorList-{colorList_id}_heatmapHeight-{heatmapHeight}_heatmapWidth-{heatmapWidth}_whatToShow-{whatToShow_id}_xAxisLabel-{xAxisLabel_id}_startLabel-{startLabel_id}_endLabel-{endLabel_id}/{filler}.pdf",
        bed="out/deepTools/plotHeatmap_sortRegions-{sortRegions}_sortUsing-{sortUsing}_averageTypeSummaryPlot-{averageTypeSummaryPlot}_colorList-{colorList_id}_heatmapHeight-{heatmapHeight}_heatmapWidth-{heatmapWidth}_whatToShow-{whatToShow_id}_xAxisLabel-{xAxisLabel_id}_startLabel-{startLabel_id}_endLabel-{endLabel_id}/{filler}.bed"
    threads:
        1
    priority:
        2
    params:
        colorList=params_deepTools_plotHeatmap_colorList,
        whatToShow=params_deepTools_plotHeatmap_whatToShow,
        xAxisLabel=params_deepTools_plotHeatmap_xAxisLabel,
        startLabel=params_deepTools_plotHeatmap_startLabel,
        endLabel=params_deepTools_plotHeatmap_endLabel,
    wildcard_constraints:
        sortRegions="descend|ascend|no",
        sortUsing="mean|median|max|min|sum|region_length",
        averageTypeSummaryPlot="mean|median|min|max|sum|std",
        colorList_id="[a-zA-Z]+",
        heatmapHeight="[0-9]+", #Plot height in cm. The default value is 28. The minimum value is 3 and the maximum is 100.
        heatmapWidth="[0-9]+", #Plot width in cm. The default value is 4 The minimum value is 1 and the maximum is 100.
        whatToShow_id="phc|ph|h|hc", # write a function to map these to path compliant ids.
        xAxisLabel_id="[a-zA-Z-]+", # write a function to map these to path compliant ids.
        startLabel_id="[a-zA-Z0-9-]+", #  write a function to map these to path compliant ids.
        endLabel_id="[a-zA-Z0-9-]+", #  write a function to map these to path compliant ids.
    shell:
        """
        {input.plotHeatmap} \
            --matrixFile {input.matrix} \
            --outFileName {output.pdf} \
            --outFileSortedRegions {output.bed} \
            --sortRegions {wildcards.sortRegions} \
            --sortUsing {wildcards.sortUsing} \
            --averageTypeSummaryPlot {wildcards.averageTypeSummaryPlot} \
            --colorList {params.colorList} \
            --heatmapHeight {wildcards.heatmapHeight} \
            --heatmapWidth {wildcards.heatmapWidth} \
            --whatToShow {params.whatToShow} \
            --xAxisLabel {params.xAxisLabel} \
            --startLabel {params.startLabel} \
            --endLabel {params.endLabel}
        """

rule deepTools_plotHeatmap_sortRegions_sortUsing_averageTypeSummaryPlot_colorList_heatmapHeight_heatmapWidth_whatToShow_xAxisLabel_startLabel_endLabel_perGroup:
    """
    Created:
        2017-10-24 10:46:31
    Test:
    """
    input:
        matrix="out/{filler}.txt.gz",
        plotHeatmap="opt/miniconda/envs/deeptools/bin/plotHeatmap"
    output:
        pdf="out/deepTools/plotHeatmap_sortRegions-{sortRegions}_sortUsing-{sortUsing}_averageTypeSummaryPlot-{averageTypeSummaryPlot}_colorList-{colorList_id}_heatmapHeight-{heatmapHeight}_heatmapWidth-{heatmapWidth}_whatToShow-{whatToShow_id}_xAxisLabel-{xAxisLabel_id}_startLabel-{startLabel_id}_endLabel-{endLabel_id}_perGroup/{filler}.pdf",
        bed="out/deepTools/plotHeatmap_sortRegions-{sortRegions}_sortUsing-{sortUsing}_averageTypeSummaryPlot-{averageTypeSummaryPlot}_colorList-{colorList_id}_heatmapHeight-{heatmapHeight}_heatmapWidth-{heatmapWidth}_whatToShow-{whatToShow_id}_xAxisLabel-{xAxisLabel_id}_startLabel-{startLabel_id}_endLabel-{endLabel_id}_perGroup/{filler}.bed"
    threads:
        1
    priority:
        2
    params:
        colorList=params_deepTools_plotHeatmap_colorList,
        whatToShow=params_deepTools_plotHeatmap_whatToShow,
        xAxisLabel=params_deepTools_plotHeatmap_xAxisLabel,
        startLabel=params_deepTools_plotHeatmap_startLabel,
        endLabel=params_deepTools_plotHeatmap_endLabel,
    wildcard_constraints:
        sortRegions="descend|ascend|no",
        sortUsing="mean|median|max|min|sum|region_length",
        averageTypeSummaryPlot="mean|median|min|max|sum|std",
        colorList_id="[a-zA-Z]+",
        heatmapHeight="[0-9]+", #Plot height in cm. The default value is 28. The minimum value is 3 and the maximum is 100.
        heatmapWidth="[0-9]+", #Plot width in cm. The default value is 4 The minimum value is 1 and the maximum is 100.
        whatToShow_id="phc|ph|h|hc", # write a function to map these to path compliant ids.
        xAxisLabel_id="[a-zA-Z-]+", # write a function to map these to path compliant ids.
        startLabel_id="[a-zA-Z0-9-]+", #  write a function to map these to path compliant ids.
        endLabel_id="[a-zA-Z0-9-]+", #  write a function to map these to path compliant ids.
    shell:
        """
        {input.plotHeatmap} \
            --matrixFile {input.matrix} \
            --outFileName {output.pdf} \
            --outFileSortedRegions {output.bed} \
            --sortRegions {wildcards.sortRegions} \
            --sortUsing {wildcards.sortUsing} \
            --averageTypeSummaryPlot {wildcards.averageTypeSummaryPlot} \
            --colorList {params.colorList} \
            --heatmapHeight {wildcards.heatmapHeight} \
            --heatmapWidth {wildcards.heatmapWidth} \
            --whatToShow {params.whatToShow} \
            --xAxisLabel {params.xAxisLabel} \
            --startLabel {params.startLabel} \
            --endLabel {params.endLabel} \
            --perGroup
        """

rule deepTools_plotHeatmap_sortRegions_sortUsing_averageTypeSummaryPlot_colorList_zMin_zMax_heatmapHeight_heatmapWidth_whatToShow_xAxisLabel_startLabel_endLabel_yMin_yMax:
    """
    Created:
        2017-08-02 22:59:39
    Test:
        out/deepTools/plotHeatmap_sortRegions-descend_sortUsing-mean_averageTypeSummaryPlot-mean_colorList-blueCyanYellowOrangeRed_zMin-0_zMax-200_heatmapHeight-28_heatmapWidth-4_whatToShow-phc_xAxisLabel-peak-center_startLabel-start_endLabel-end_yMin-0_yMax-100/deepTools/computeMatrix_scale-regions_beforeRegionStartLength-5000_afterRegionStartLength-5000_bed-mm10-h4k5ac-bu-peaks-gfold-diff-signif_bw-mm10-MNase-wt-h4k5k8acbu.pdf
    """
    input:
        matrix="out/{filler}.txt.gz",
        plotHeatmap="opt/miniconda/envs/deeptools/bin/plotHeatmap"
    output:
        pdf="out/deepTools/plotHeatmap_sortRegions-{sortRegions}_sortUsing-{sortUsing}_averageTypeSummaryPlot-{averageTypeSummaryPlot}_colorList-{colorList_id}_zMin-{zMin}_zMax-{zMax}_heatmapHeight-{heatmapHeight}_heatmapWidth-{heatmapWidth}_whatToShow-{whatToShow_id}_xAxisLabel-{xAxisLabel_id}_startLabel-{startLabel_id}_endLabel-{endLabel_id}_yMin-{yMin}_yMax-{yMax}/{filler}.pdf",
        bed="out/deepTools/plotHeatmap_sortRegions-{sortRegions}_sortUsing-{sortUsing}_averageTypeSummaryPlot-{averageTypeSummaryPlot}_colorList-{colorList_id}_zMin-{zMin}_zMax-{zMax}_heatmapHeight-{heatmapHeight}_heatmapWidth-{heatmapWidth}_whatToShow-{whatToShow_id}_xAxisLabel-{xAxisLabel_id}_startLabel-{startLabel_id}_endLabel-{endLabel_id}_yMin-{yMin}_yMax-{yMax}/{filler}.bed"
    threads:
        1
    priority:
        2
    params:
        colorList=params_deepTools_plotHeatmap_colorList,
        whatToShow=params_deepTools_plotHeatmap_whatToShow,
        xAxisLabel=params_deepTools_plotHeatmap_xAxisLabel,
        startLabel=params_deepTools_plotHeatmap_startLabel,
        endLabel=params_deepTools_plotHeatmap_endLabel,
    wildcard_constraints:
        sortRegions="descend|ascend|no",
        sortUsing="mean|median|max|min|sum|region_length",
        averageTypeSummaryPlot="mean|median|min|max|sum|std",
        colorList_id="[a-zA-Z]+",
        zMin="[0-9-\.]+",
        zMax="[0-9-\.]+",
        heatmapHeight="[0-9]+", #Plot height in cm. The default value is 28. The minimum value is 3 and the maximum is 100.
        heatmapWidth="[0-9]+", #Plot width in cm. The default value is 4 The minimum value is 1 and the maximum is 100.
        whatToShow_id="phc|ph|h|hc", # write a function to map these to path compliant ids.
        xAxisLabel_id="[a-zA-Z-]+", # write a function to map these to path compliant ids.
        startLabel_id="[a-zA-Z0-9-]+", #  write a function to map these to path compliant ids.
        endLabel_id="[a-zA-Z0-9-]+", #  write a function to map these to path compliant ids.
        yMin="[0-9-\.]+",
        yMax="[0-9-\.]+",
    shell:
        """
        {input.plotHeatmap} \
            --matrixFile {input.matrix} \
            --outFileName {output.pdf} \
            --outFileSortedRegions {output.bed} \
            --sortRegions {wildcards.sortRegions} \
            --sortUsing {wildcards.sortUsing} \
            --averageTypeSummaryPlot {wildcards.averageTypeSummaryPlot} \
            --colorList {params.colorList} \
            --zMin '{wildcards.zMin}' \
            --zMax '{wildcards.zMax}' \
            --heatmapHeight {wildcards.heatmapHeight} \
            --heatmapWidth {wildcards.heatmapWidth} \
            --whatToShow {params.whatToShow} \
            --xAxisLabel {params.xAxisLabel} \
            --startLabel {params.startLabel} \
            --endLabel {params.endLabel} \
            --yMin '{wildcards.yMin}' \
            --yMax '{wildcards.yMax}'
        """

rule deepTools_plotHeatmap_kmeans_sortRegions_sortUsing_averageTypeSummaryPlot_colorList_heatmapHeight_heatmapWidth_whatToShow_xAxisLabel_startLabel_endLabel:
    """
    Created:
        2017-10-19 10:05:10
    Test:
    """
    input:
        matrix="out/{filler}.txt.gz",
        plotHeatmap="opt/miniconda/envs/deeptools/bin/plotHeatmap"
    output:
        pdf="out/deepTools/plotHeatmap_kmeans-{kmeans}_sortRegions-{sortRegions}_sortUsing-{sortUsing}_averageTypeSummaryPlot-{averageTypeSummaryPlot}_colorList-{colorList_id}_heatmapHeight-{heatmapHeight}_heatmapWidth-{heatmapWidth}_whatToShow-{whatToShow_id}_xAxisLabel-{xAxisLabel_id}_startLabel-{startLabel_id}_endLabel-{endLabel_id}/{filler}.pdf",
        bed="out/deepTools/plotHeatmap_kmeans-{kmeans}_sortRegions-{sortRegions}_sortUsing-{sortUsing}_averageTypeSummaryPlot-{averageTypeSummaryPlot}_colorList-{colorList_id}_heatmapHeight-{heatmapHeight}_heatmapWidth-{heatmapWidth}_whatToShow-{whatToShow_id}_xAxisLabel-{xAxisLabel_id}_startLabel-{startLabel_id}_endLabel-{endLabel_id}/{filler}.bed"
    threads:
        1
    priority:
        2
    params:
        colorList=params_deepTools_plotHeatmap_colorList,
        whatToShow=params_deepTools_plotHeatmap_whatToShow,
        xAxisLabel=params_deepTools_plotHeatmap_xAxisLabel,
        startLabel=params_deepTools_plotHeatmap_startLabel,
        endLabel=params_deepTools_plotHeatmap_endLabel,
    wildcard_constraints:
        kmeans="[0-9]+",  #To check: kclasses=1 should be equivalent to 'None' and avoid the if/else below.
        sortRegions="descend|ascend|no",
        sortUsing="mean|median|max|min|sum|region_length",
        averageTypeSummaryPlot="mean|median|min|max|sum|std",
        colorList_id="[a-zA-Z]+",
        heatmapHeight="[0-9]+", #Plot height in cm. The default value is 28. The minimum value is 3 and the maximum is 100.
        heatmapWidth="[0-9]+", #Plot width in cm. The default value is 4 The minimum value is 1 and the maximum is 100.
        whatToShow_id="phc|ph|h|hc", # write a function to map these to path compliant ids.
        xAxisLabel_id="[a-zA-Z-]+", # write a function to map these to path compliant ids.
        startLabel_id="[a-zA-Z0-9-]+", #  write a function to map these to path compliant ids.
        endLabel_id="[a-zA-Z0-9-]+", #  write a function to map these to path compliant ids.
    shell:
        """
        {input.plotHeatmap} \
            --matrixFile {input.matrix} \
            --outFileName {output.pdf} \
            --outFileSortedRegions {output.bed} \
            --kmeans {wildcards.kmeans} \
            --sortRegions {wildcards.sortRegions} \
            --sortUsing {wildcards.sortUsing} \
            --averageTypeSummaryPlot {wildcards.averageTypeSummaryPlot} \
            --colorList {params.colorList} \
            --heatmapHeight {wildcards.heatmapHeight} \
            --heatmapWidth {wildcards.heatmapWidth} \
            --whatToShow {params.whatToShow} \
            --xAxisLabel {params.xAxisLabel} \
            --startLabel {params.startLabel} \
            --endLabel {params.endLabel}
        """

rule deepTools_plotHeatmap_kmeans_sortRegions_sortUsing_averageTypeSummaryPlot_colorList_zMin_zMax_heatmapHeight_heatmapWidth_whatToShow_xAxisLabel_startLabel_endLabel_yMin_yMax:
    """
    Created:
        2017-05-07 13:19:39
    Test:
        out/deepTools/plotHeatmap_kmeans-20_sortRegions-descend_sortUsing-mean_averageTypeSummaryPlot-mean_colorList-blueCyanYellowOrangeRed_zMin-0_zMax-200_heatmapHeight-28_heatmapWidth-4_whatToShow-phc_xAxisLabel-peak-center_startLabel-start_endLabel-end_yMin-0_yMax-100/deepTools/computeMatrix_scale-regions_beforeRegionStartLength-5000_afterRegionStartLength-5000_bed-mm10-h2alap1-peaks_bw-mm10-MNase-wt.pdf
    """
    input:
        matrix="out/{filler}.txt.gz",
        plotHeatmap="opt/miniconda/envs/deeptools/bin/plotHeatmap"
    output:
        pdf="out/deepTools/plotHeatmap_kmeans-{kmeans}_sortRegions-{sortRegions}_sortUsing-{sortUsing}_averageTypeSummaryPlot-{averageTypeSummaryPlot}_colorList-{colorList_id}_zMin-{zMin}_zMax-{zMax}_heatmapHeight-{heatmapHeight}_heatmapWidth-{heatmapWidth}_whatToShow-{whatToShow_id}_xAxisLabel-{xAxisLabel_id}_startLabel-{startLabel_id}_endLabel-{endLabel_id}_yMin-{yMin}_yMax-{yMax}/{filler}.pdf",
        bed="out/deepTools/plotHeatmap_kmeans-{kmeans}_sortRegions-{sortRegions}_sortUsing-{sortUsing}_averageTypeSummaryPlot-{averageTypeSummaryPlot}_colorList-{colorList_id}_zMin-{zMin}_zMax-{zMax}_heatmapHeight-{heatmapHeight}_heatmapWidth-{heatmapWidth}_whatToShow-{whatToShow_id}_xAxisLabel-{xAxisLabel_id}_startLabel-{startLabel_id}_endLabel-{endLabel_id}_yMin-{yMin}_yMax-{yMax}/{filler}.bed"
    threads:
        1
    priority:
        2
    params:
        colorList=params_deepTools_plotHeatmap_colorList,
        whatToShow=params_deepTools_plotHeatmap_whatToShow,
        xAxisLabel=params_deepTools_plotHeatmap_xAxisLabel,
        startLabel=params_deepTools_plotHeatmap_startLabel,
        endLabel=params_deepTools_plotHeatmap_endLabel,
    wildcard_constraints:
        kmeans="[0-9]+",  #To check: kclasses=1 should be equivalent to 'None' and avoid the if/else below.
        sortRegions="descend|ascend|no",
        sortUsing="mean|median|max|min|sum|region_length",
        averageTypeSummaryPlot="mean|median|min|max|sum|std",
        colorList_id="[a-zA-Z]+",
        zMin="[0-9-\.]+",
        zMax="[0-9-\.]+",
        heatmapHeight="[0-9]+", #Plot height in cm. The default value is 28. The minimum value is 3 and the maximum is 100.
        heatmapWidth="[0-9]+", #Plot width in cm. The default value is 4 The minimum value is 1 and the maximum is 100.
        whatToShow_id="phc|ph|h|hc", # write a function to map these to path compliant ids.
        xAxisLabel_id="[a-zA-Z-]+", # write a function to map these to path compliant ids.
        startLabel_id="[a-zA-Z0-9-]+", #  write a function to map these to path compliant ids.
        endLabel_id="[a-zA-Z0-9-]+", #  write a function to map these to path compliant ids.
        yMin="[0-9-\.]+",
        yMax="[0-9-\.]+",
    shell:
        """
        {input.plotHeatmap} \
            --matrixFile {input.matrix} \
            --outFileName {output.pdf} \
            --outFileSortedRegions {output.bed} \
            --kmeans {wildcards.kmeans} \
            --sortRegions {wildcards.sortRegions} \
            --sortUsing {wildcards.sortUsing} \
            --averageTypeSummaryPlot {wildcards.averageTypeSummaryPlot} \
            --colorList {params.colorList} \
            --zMin '{wildcards.zMin}' \
            --zMax '{wildcards.zMax}' \
            --heatmapHeight {wildcards.heatmapHeight} \
            --heatmapWidth {wildcards.heatmapWidth} \
            --whatToShow {params.whatToShow} \
            --xAxisLabel {params.xAxisLabel} \
            --startLabel {params.startLabel} \
            --endLabel {params.endLabel} \
            --yMin '{wildcards.yMin}' \
            --yMax '{wildcards.yMax}'
        """

rule deepTools_plotHeatmap_kmeans_sortRegions_sortUsing_averageTypeSummaryPlot_colorList_zMin_zMax_heatmapHeight_heatmapWidth_whatToShow_xAxisLabel_refPointLabel_yMin_yMax:
    """
    Created:
        2017-05-07 13:19:39
    Test:
        out/deepTools/plotHeatmap_kmeans-1_sortRegions-descend_sortUsing-mean_averageTypeSummaryPlot-mean_colorList-blueCyanYellowOrangeRed_zMin-0_zMax-50_heatmapHeight-28_heatmapWidth-4_whatToShow-phc_xAxisLabel-distance-to-nucleosome-center-bp_refPointLabel-0_yMin-0_yMax-50/deepTools/computeMatrix_reference-point_referencePoint-center_beforeRegionStartLength-100_afterRegionStartLength-100_bed-mm9-nuc-c-wt-danpos-maxfuz-40-minsmt-200_bw-MNase-wt-2016-12-02.pdf
        out/deepTools/plotHeatmap_kmeans-3_sortRegions-descend_sortUsing-mean_averageTypeSummaryPlot-mean_colorList-blueCyanYellowOrangeRed_zMin-0_zMax-50_heatmapHeight-28_heatmapWidth-4_whatToShow-phc_xAxisLabel-distance-to-nucleosome-center-bp_refPointLabel-0_yMin-0_yMax-50/deepTools/computeMatrix_reference-point_referencePoint-center_beforeRegionStartLength-200_afterRegionStartLength-200_bed-mm9-nuc-c-wt-danpos-maxfuz-50-minsmt-400_bw-MNase-wt-2017-05-08.pdf
    """
    input:
        matrix="out/{filler}.txt.gz",
        plotHeatmap="opt/miniconda/envs/deeptools/bin/plotHeatmap"
    output:
        pdf="out/deepTools/plotHeatmap_kmeans-{kmeans}_sortRegions-{sortRegions}_sortUsing-{sortUsing}_averageTypeSummaryPlot-{averageTypeSummaryPlot}_colorList-{colorList_id}_zMin-{zMin}_zMax-{zMax}_heatmapHeight-{heatmapHeight}_heatmapWidth-{heatmapWidth}_whatToShow-{whatToShow_id}_xAxisLabel-{xAxisLabel_id}_refPointLabel-{refPointLabel_id}_yMin-{yMin}_yMax-{yMax}/{filler}.pdf",
        bed="out/deepTools/plotHeatmap_kmeans-{kmeans}_sortRegions-{sortRegions}_sortUsing-{sortUsing}_averageTypeSummaryPlot-{averageTypeSummaryPlot}_colorList-{colorList_id}_zMin-{zMin}_zMax-{zMax}_heatmapHeight-{heatmapHeight}_heatmapWidth-{heatmapWidth}_whatToShow-{whatToShow_id}_xAxisLabel-{xAxisLabel_id}_refPointLabel-{refPointLabel_id}_yMin-{yMin}_yMax-{yMax}/{filler}.bed"
    threads:
        1
    priority:
        2
    params:
        colorList=params_deepTools_plotHeatmap_colorList,
        whatToShow=params_deepTools_plotHeatmap_whatToShow,
        xAxisLabel=params_deepTools_plotHeatmap_xAxisLabel,
        refPointLabel=params_deepTools_plotHeatmap_refPointLabel,
    wildcard_constraints:
        kmeans="[0-9]+",  #To check: kclasses=1 should be equivalent to 'None' and avoid the if/else below.
        sortRegions="descend|ascend|no",
        sortUsing="mean|median|max|min|sum|region_length",
        averageTypeSummaryPlot="mean|median|min|max|sum|std",
        colorList_id="[a-zA-Z]+",
        zMin="[0-9-\.]+",
        zMax="[0-9-\.]+",
        heatmapHeight="[0-9]+", #Plot height in cm. The default value is 28. The minimum value is 3 and the maximum is 100.
        heatmapWidth="[0-9]+", #Plot width in cm. The default value is 4 The minimum value is 1 and the maximum is 100.
        whatToShow_id="phc|ph|h|hc", # write a function to map these to path compliant ids.
        xAxisLabel_id="[a-zA-Z-]+", # write a function to map these to path compliant ids.
        refPointLabel_id="[a-zA-Z0-9-]+", #  write a function to map these to path compliant ids.
        yMin="[0-9-\.]+",
        yMax="[0-9-\.]+",
    shell:
        """
        {input.plotHeatmap} \
            --matrixFile {input.matrix} \
            --outFileName {output.pdf} \
            --outFileSortedRegions {output.bed} \
            --kmeans {wildcards.kmeans} \
            --sortRegions {wildcards.sortRegions} \
            --sortUsing {wildcards.sortUsing} \
            --averageTypeSummaryPlot {wildcards.averageTypeSummaryPlot} \
            --colorList {params.colorList} \
            --zMin '{wildcards.zMin}' \
            --zMax '{wildcards.zMax}' \
            --heatmapHeight {wildcards.heatmapHeight} \
            --heatmapWidth {wildcards.heatmapWidth} \
            --whatToShow {params.whatToShow} \
            --xAxisLabel {params.xAxisLabel} \
            --refPointLabel {params.refPointLabel} \
            --yMin '{wildcards.yMin}' \
            --yMax '{wildcards.yMax}'
        """

rule deepTools_plotHeatmap_hclust_sortRegions_sortUsing_averageTypeSummaryPlot_colorList_zMin_zMax_heatmapHeight_heatmapWidth_whatToShow_xAxisLabel_refPointLabel_yMin_yMax:
    """
    Created:
        2017-05-08 16:29:27
    Test:
        out/deepTools/plotHeatmap_hclust-5_sortRegions-descend_sortUsing-region_length_averageTypeSummaryPlot-mean_colorList-blueCyanYellowOrangeRed_zMin-0_zMax-200_heatmapHeight-28_heatmapWidth-4_whatToShow-phc_xAxisLabel-peak-center_refPointLabel-0_yMin-0_yMax-100/deepTools/computeMatrix_reference-point_referencePoint-center_beforeRegionStartLength-1000_afterRegionStartLength-1000_bed-mm10-h2alap1-peaks_bw-mm10-MNase-wt.pdf
    """
    input:
        matrix="out/{filler}.txt.gz",
        plotHeatmap="opt/miniconda/envs/deeptools/bin/plotHeatmap"
    output:
        pdf="out/deepTools/plotHeatmap_hclust-{hclust}_sortRegions-{sortRegions}_sortUsing-{sortUsing}_averageTypeSummaryPlot-{averageTypeSummaryPlot}_colorList-{colorList_id}_zMin-{zMin}_zMax-{zMax}_heatmapHeight-{heatmapHeight}_heatmapWidth-{heatmapWidth}_whatToShow-{whatToShow_id}_xAxisLabel-{xAxisLabel_id}_refPointLabel-{refPointLabel_id}_yMin-{yMin}_yMax-{yMax}/{filler}.pdf",
        bed="out/deepTools/plotHeatmap_hclust-{hclust}_sortRegions-{sortRegions}_sortUsing-{sortUsing}_averageTypeSummaryPlot-{averageTypeSummaryPlot}_colorList-{colorList_id}_zMin-{zMin}_zMax-{zMax}_heatmapHeight-{heatmapHeight}_heatmapWidth-{heatmapWidth}_whatToShow-{whatToShow_id}_xAxisLabel-{xAxisLabel_id}_refPointLabel-{refPointLabel_id}_yMin-{yMin}_yMax-{yMax}/{filler}.bed"
    threads:
        1
    priority:
        2
    params:
        colorList=params_deepTools_plotHeatmap_colorList,
        whatToShow=params_deepTools_plotHeatmap_whatToShow,
        xAxisLabel=params_deepTools_plotHeatmap_xAxisLabel,
        refPointLabel=params_deepTools_plotHeatmap_refPointLabel,
    wildcard_constraints:
        hclust="[0-9]+",  #To check: kclasses=1 should be equivalent to 'None' and avoid the if/else below.
        sortRegions="descend|ascend|no",
        sortUsing="mean|median|max|min|sum|region_length",
        averageTypeSummaryPlot="mean|median|min|max|sum|std",
        colorList_id="[a-zA-Z]+",
        zMin="[0-9-\.]+",
        zMax="[0-9-\.]+",
        heatmapHeight="[0-9]+", #Plot height in cm. The default value is 28. The minimum value is 3 and the maximum is 100.
        heatmapWidth="[0-9]+", #Plot width in cm. The default value is 4 The minimum value is 1 and the maximum is 100.
        whatToShow_id="phc|ph|h|hc", # write a function to map these to path compliant ids.
        xAxisLabel_id="[a-zA-Z-]+", # write a function to map these to path compliant ids.
        refPointLabel_id="[a-zA-Z0-9-]+", #  write a function to map these to path compliant ids.
        yMin="[0-9-\.]+",
        yMax="[0-9-\.]+",
    shell:
        """
        {input.plotHeatmap} \
            --matrixFile {input.matrix} \
            --outFileName {output.pdf} \
            --outFileSortedRegions {output.bed} \
            --hclust {wildcards.hclust} \
            --sortRegions {wildcards.sortRegions} \
            --sortUsing {wildcards.sortUsing} \
            --averageTypeSummaryPlot {wildcards.averageTypeSummaryPlot} \
            --colorList {params.colorList} \
            --zMin '{wildcards.zMin}' \
            --zMax '{wildcards.zMax}' \
            --heatmapHeight {wildcards.heatmapHeight} \
            --heatmapWidth {wildcards.heatmapWidth} \
            --whatToShow {params.whatToShow} \
            --xAxisLabel {params.xAxisLabel} \
            --refPointLabel {params.refPointLabel} \
            --yMin '{wildcards.yMin}' \
            --yMax '{wildcards.yMax}'
        """
