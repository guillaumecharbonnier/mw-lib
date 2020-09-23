# Including global constraints on wildcards
# Including constraints allow easier developping, debugging and decrease DAG resolution time by reducing the complexity of combinations between wildcards.

# "[\w_-]+" short for "[a-zA-Z0-9_-]+"

wildcard_constraints:
    extra = "[^\/]*",
    #extra = "[a-zA-Z0-9-_:.]*", # extra should be allowed to be empty. (. and : for ffmpeg)
    #extra = ".*",
    # IDs
    # By convention, all wildcards ending with '_id' should be used as key for python functions.
    # To prevent the hassle of defining them individually like for chrominfo_id and fa_genome_id, I have decided they have to follow these constraints below:
    bed_id = "bed-[a-zA-Z0-9-]+",
    b_id = "[a-zA-Z0-9-]+",
    bed_list_id = "bed-[a-zA-Z0-9-]+",
    bw_id  = "bw-[a-zA-Z0-9-]+",
    gff_id = "gff-[a-zA-Z0-9-]+",
    gtf_id = "gtf-[a-zA-Z0-9-]+",
    vcfgz_id = "vcfgz-[a-zA-Z0-9-]+",
    vcf_id = "vcf-[a-zA-Z0-9-]+",
    ebwt_id = "ebwt-[a-zA-Z0-9-]+",
    ref_id = "[a-zA-Z0-9-]+",
    chrominfo_id  = "chrominfo-[a-zA-Z0-9-]+",
    fa_genome_id = "fa-[a-zA-Z0-9-]+",
    bwa_index_id = "bwa-index-[a-zA-Z0-9-]+",
    kallisto_idx_id = "kallisto-idx-[a-zA-Z0-9-]+",
    chain_id = "chain-[a-zA-Z0-9-]+",
    basename = "[a-zA-Z0-9-_]+",
    filename = "[a-zA-Z0-9_-]+",
    #index = "mm9|mm10|hg18|hg19|hg38|GRCh38|GRCm38|dm6|BDGP6|sacCer3|GRCh38-ensembl-r100", # Reference genome assembly
    index = "[a-zA-Z0-9-]+",
    assembly = "mm9|mm10|hg18|hg19|hg38|GRCh38|GRCm38|dm6|BDGP6|sacCer3|GRCh38-ensembl-r100", # Reference genome assembly
    #gff_id = "mm9|mm10|mm10-chr1|hg18|hg19|hg38|GRCh38|GRCm38|dm6|sacCer3|GRCm38-Ensembl-91-chr19", # Reference genome assembly
    runtype = "se|pe", # Used to discriminate single and paired-ends specific rules.
    refType = "TSS|TES|center", # wildcards for deepTools computeMatrix, reference used in reference-point mode.
    signal = "bamCoverage|bamCompare", #wildcards for bw signal used by deepTools computeMatrix.
    lengthAround = "[0-9]+", #wildcard for deepTools computeMatrix, to define the region evaluated around the reference.
    geo_experiment = "GSE[0-9]+",
    sra_sample = "SRR[0-9]+",
    se_pe="se|pe",
    mate_prefix = "R|", # Sometimes mates are identified by 'R1/R2' or just '1/2'...
    seed = "[0-9]+",
    seed1 = "[0-9]+",
    seed2 = "[0-9]+",
    ppr = "positions|regions|peaks",
    pseudocount = "0|1|0.000001",
    #legacy below:
    sickle_qual_threshold = "[0-9]+",
    sickle_qual_type = "illumina|sanger|solexa",
    # deepTools plotProfile and plotHeatmap
    deepTools_extendReads = "|[0-9]+",
    deepTools_plotCorrelation_corMethod = "spearman|pearson",
    deepTools_plotCorrelation_whatToPlot = "heatmap|scatterplot",
    deepTools_computeMatrix_referencePoint = "TSS|TES|center",
    deepTools_computeMatrix_beforeRegionStartLength = "[0-9]+",
    deepTools_computeMatrix_afterRegionStartLength = "[0-9]+",
    deepTools_computeMatrix_binSize = "[0-9]+",
    deepTools_plotHeatmap_kmeans = "[0-9]+",
    deepTools_plotHeatmap_sortRegions = "descend|ascend|no",
    deepTools_plotHeatmap_sortUsing = "mean|median|max|min|sum|region_length",
    deepTools_plotHeatmap_averageTypeSummaryPlot = "mean|median|min|max|sum|std",
    deepTools_plotHeatmap_missingDataColor = "0|1",
    deepTools_plotHeatmap_colorList_id = "[a-zA-Z]+",
    deepTools_plotHeatmap_zMin = "[0-9]+",
    deepTools_plotHeatmap_zMax = "[0-9]+",
    deepTools_plotHeatmap_heatmapHeight = "[0-9]+", #Plot height in cm. The default value is 28. The minimum value is 3 and the maximum is 100.
    deepTools_plotHeatmap_heatmapWidth = "[0-9]+", #Plot width in cm. The default value is 4 The minimum value is 1 and the maximum is 100.
    deepTools_plotHeatmap_whatToShow_id = "phc|ph|h|hc", # write a function to map these to path compliant ids.
    deepTools_plotHeatmap_xAxisLabel_id = "[a-zA-Z-]+", # write a function to map these to path compliant ids.
    deepTools_plotHeatmap_refPointLabel_id = "[a-zA-Z0-9-]+", #  write a function to map these to path compliant ids.
    deepTools_plotHeatmap_boxAroundHeatmaps = "yes|no"
    # CapStarr
