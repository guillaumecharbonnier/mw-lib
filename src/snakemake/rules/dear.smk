DEAR_OUTFILES=[
"DESeq2_CONvsTRT_DEG.txt",
"DESeq2_CONvsTRT_differential_expression.txt",
"DESeq2_CONvsTRT_Histogram_pvalues.jpeg",
"DESeq2_CONvsTRT_Volcanoplot.jpeg",
"DESeq2_Dispersionplot.jpeg",
"DESeq2_genesconsidered.txt",
"DESeq2_MAplot.jpeg",
"DESeq2_PCA_scree.jpeg",
"DESeq2_SamplesHeatMap_normalizedcounts.jpeg",
"DESeq2_vst_normalizedcounts.txt",
"edgeR_DEG.txt",
"edgeR_differential_expression.txt",
"edgeR_dispersionPlot.jpeg",
"edgeR_genesconsidered.txt",
"edgeR_MAplot.jpeg",
"edgeR_normalizedCounts_cpm.txt",
"edgeR_PCA_scree.jpeg",
"edgeR_Volcanoplot.jpeg",
"FeatureCounts_CountingStatistics.txt",
"FeatureCounts.log",
"FeatureCounts_NoFeaturePlot.jpeg",
"FeatureCounts_RawCounts.RData",
"FeatureCounts_RawCounts.txt",
"Limma-voom_consideredgenes.txt",
"Limma-voom_DEG.txt",
"Limma-voom_differential_expression.txt",
"Limma-voom_normalizedcounts.txt",
"Limma-voom_PCA_scree.jpeg",
"Limma-voom_Volcanoplot.jpeg",
"VennDiagram-DEG.jpeg"]

rule dear:
    """
    Created:
        2019-02-05 22:28:21
    Aim:
    Todo:
        Bam could be replaced by a function that parse the TSV looking for the bam.
    Doc:
        https://github.com/wdecoster/DEA.R
    Test:
        out/dear/dear-hg38-casero2016-thy3-thy4_gtf-hg38-ensembl_bam-hg38-casero2016-thy3-thy4/done
    """
    input:
        dear="../DEA.R/DEA/DEA.R",
        bam = lambda wildcards: eval(config['ids'][wildcards.bam_list_id]),
        tsv = lambda wildcards: config['ids'][wildcards.dear_id],
        gtf = lambda wildcards: config['ids'][wildcards.gtf_id]
    output:
        expand("out/dear/{{dear_id}}_{{gtf_id}}_{{bam_list_id}}/{file}", file=DEAR_OUTFILES)
    params:
        outdir="out/dear/{dear_id}_{gtf_id}_{bam_list_id}"
    conda:
        "../envs/dear.yaml"
    shell:
        "{input.dear} {input.tsv} {input.gtf} --output {params.outdir}"
