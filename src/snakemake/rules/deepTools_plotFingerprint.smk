rule deepTools_plotFingerprint_extra_lambda:
    """
    Created:
        2018-01-18 15:08:41
    Aim:
        QC for ChIP-Seq. This "lambda" version of the rule allows to plot multiple samples on the same plot.
    Doc:
        https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html
        https://deeptools.readthedocs.io/en/develop/content/feature/plotFingerprint_QC_metrics.html
    Note:
        --smartLabels nearly always useful.
        --outRawCounts added to satisfy multiqc: https://multiqc.info/docs/#deeptools
    Test:
        out/deepTools/plotFingerprint_extendReads-100/Blueprint-thymic-populations-H3K27ac.pdf
        out/deepTools/plotFingerprint_extendReads-100/Blueprint-thymic-populations-H3K4me1.pdf
        out/deepTools/plotFingerprint_extendReads-/mm10_H4K5ac_comparison.pdf
        out/deepTools/plotFingerprint_bam-hg38-run246.pdf
    """
    input:
        bam = lambda wildcards: eval(mwconf['ids'][wildcards.bam_list_id]),
        bai = lambda wildcards: [path + '.bai' for path in eval(mwconf['ids'][wildcards.bam_list_id])]
    output:
        pdf = "out/{tool}{extra}_{bam_list_id}.pdf",
        txt = "out/{tool}{extra}_{bam_list_id}.txt"
    log:
              "out/{tool}{extra}_{bam_list_id}.log"
    benchmark:
              "out/{tool}{extra}_{bam_list_id}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="deepTools/plotFingerprint"
    conda:
        "../envs/deeptools.yaml"
    threads:
        MAX_THREADS
    shell:
        "plotFingerprint --bamfiles {input.bam} --plotFile {output.pdf} --outQualityMetrics {output.txt} --numberOfProcessors {threads} {params.extra} &> {log}"

rule deepTools_plotFingerprint_extra_filler:
    """
    Created:
        2018-01-18 15:08:41
    Aim:
        QC for ChIP-Seq. This "filler" version of the rule is convenient for one sample.
    Doc:
        https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html
        https://deeptools.readthedocs.io/en/develop/content/feature/plotFingerprint_QC_metrics.html
    Note:
        --smartLabels nearly always useful.
    Test:
        out/deepTools/plotFingerprint/samtools/index/samtools/sort/samtools/view_sam_to_bam_-q_30/bowtie2/se_mm10/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR3126243.metrics.tsv
    """
    input:
        bam = "out/{filler}.bam",
        bai = "out/{filler}.bam.bai"
    output:
        pdf         = "out/{tool}{extra}/{filler}.pdf",
        tsv_metrics = "out/{tool}{extra}/{filler}.metrics.tsv",
        tsv_counts  = "out/{tool}{extra}/{filler}.counts.tsv"
    log:
                      "out/{tool}{extra}/{filler}.log"
    benchmark:
                      "out/{tool}{extra}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="deepTools/plotFingerprint"
    conda:
        "../envs/deeptools.yaml"
    threads:
        MAX_THREADS
    shell:
        "plotFingerprint --bamfiles {input.bam} "
        "--plotFile {output.pdf} "
        "--outQualityMetrics {output.tsv_metrics} "
        "--outRawCounts {output.tsv_counts} "
        "--numberOfProcessors {threads} "
        "{params.extra} &> {log}"


