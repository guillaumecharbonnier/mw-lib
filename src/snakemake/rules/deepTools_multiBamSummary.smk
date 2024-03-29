# out/sort/coordinates_bed/cat/hg38-macs2-peaks-H3K27ac-thymus.bed

rule deepTools_multiBamSummary_BED_extra:
    """
    Created:
        2016-08-24 15h07
    Aim:
        Create a matrix of read coverages for genomic regions used by other deepTools functions.
    Doc:
        https://deeptools.readthedocs.io/en/develop/content/tools/multiBamSummary.html
    Test:
        #out/deepTools/multiBamSummary_BED_-e_100/bed-hg38-macs2-peaks-H3K27ac-thymus/bam-Blueprint-thymic-populations-H3K27ac.npz
    """
    input:
        bam = lambda wildcards: eval(mwconf['ids'][wildcards.bam_list_id]),
        bai = lambda wildcards: [path + '.bai' for path in eval(mwconf['ids'][wildcards.bam_list_id])],
        bed = lambda wildcards: eval(mwconf['ids'][wildcards.bed_id])
    output:
        npz="out/{tool}{extra}/{bed_id}/{bam_list_id}.npz"
    log:
            "out/{tool}{extra}/{bed_id}/{bam_list_id}.log"
    benchmark:
            "out/{tool}{extra}/{bed_id}/{bam_list_id}.benchmark.tsv"
    params:
        extra = params_extra
    threads:
        MAX_THREADS
    wildcard_constraints:
        tool="deepTools/multiBamSummary_BED"
    conda:
        "../envs/deeptools.yaml"
    shell:
        "multiBamSummary BED-file --bamfiles {input.bam} --BED {input.bed} {params.extra} --outFileName {output.npz} --numberOfProcessors {threads} &> {log}"

rule deepTools_multiBamSummary_bins_extra:
    """
    Created:
        2016-08-24 15h07
    Aim:
        Create a matrix of read coverages for genomic regions used by other deepTools functions.
    Doc:
        https://deeptools.readthedocs.io/en/develop/content/tools/multiBamSummary.html
    Test:
        out/deepTools/multiBamSummary_bins_-e_147/bam-mm10-H4K5K8acbu-R.npz
    """
    input:
        bam = lambda wildcards: eval(mwconf['ids'][wildcards.bam_list_id]),
        bai = lambda wildcards: [path + '.bai' for path in eval(mwconf['ids'][wildcards.bam_list_id])]
    output:
        npz="out/{tool}{extra}/{bam_list_id}.npz"
    log:
            "out/{tool}{extra}/{bam_list_id}.log"
    benchmark:
            "out/{tool}{extra}/{bam_list_id}.benchmark.tsv"
    params:
        extra = params_extra
    threads:
        MAX_THREADS
    wildcard_constraints:
        tool="deepTools/multiBamSummary_bins",
        bam_list_id="bam-[A-Za-z0-9-]+"
    conda:
        "../envs/deeptools.yaml"
    shell:
        "multiBamSummary bins --bamfiles {input.bam} {params.extra} --outFileName {output.npz} --numberOfProcessors {threads} &> {log}"


rule deepTools_plotPCA:
    """
    Created:
        2018-01-18 15:23:47
    Doc:
        https://deeptools.readthedocs.io/en/develop/content/tools/plotPCA.html
    Note:
        --plotWidth 20 --plotHeight 20\
    Test:
        out/deepTools/plotPCA/deepTools/multiBamSummary_bins_extendReads-100/Blueprint-thymic-populations-H3K27ac.pdf
        out/deepTools/plotPCA/deepTools/multiBamSummary_BED-hg38-macs2-peaks-H3K27ac-thymus_extendReads-100/Blueprint-thymic-populations-H3K27ac.pdf
    """
    input:
        npz="out/{filler}.npz"
    output:
        pdf="out/{tool}{extra}/{filler}.pdf"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="deepTools/plotPCA"
    conda:
        "../envs/deeptools.yaml"
    shell:
        "plotPCA --corData {input.npz} --plotFile {output.pdf} {params.extra}"

rule deepTools_plotCorrelation_extra:
    """
    Created:
        2016-08-24 15h44
    Doc:
        https://deeptools.readthedocs.io/en/develop/content/tools/plotCorrelation.html
    Test:
        out/deepTools/plotCorrelation_corMethod-spearman_whatToPlot-scatterplot/deepTools/multiBamSummary_bins_extendReads-100/Blueprint-thymic-populations-H3K27ac.pdf
        ChIP-Seq_Spike-in_REH_pe_GRCh38_--very-sensitive-local
    """
    input:
        npz="out/{filler}.npz"
    output: 
        pdf="out/{tool}{extra}_corMethod-{corMethod}_whatToPlot-{whatToPlot}/{filler}.{ext}"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="deepTools/plotCorrelation",
        corMethod="spearman|pearson",
        whatToPlot="heatmap|scatterplot",
        ext="pdf|png|svg|eps"
    shell:
        "plotCorrelation --corData {input.npz} --plotFile {output.pdf} --corMethod {wildcards.corMethod} --whatToPlot {wildcards.whatToPlot} {params.extra}"

