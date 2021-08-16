rule htseq_count_extra:
    """
    Aim:
        Counts
    Doc:
        https://htseq.readthedocs.io/en/master/count.html
    Tests:
        out/htseq/count_-r_name_-s_no_--nonunique_all_gtf-GRCh38-ensembl_bam-hg38-Casero2016-thy3-thy4.tsv
    """
    input:
        gtf_saf = lambda wildcards: eval(mwconf['ids'][wildcards.gtf_saf_id]),
        bam     = lambda wildcards: eval(mwconf['ids'][wildcards.bam_list_id])
    output:
        tsv = "out/{tool}{extra}_{gtf_saf_id}_{bam_list_id}.tsv",
        txt = "out/{tool}{extra}_{gtf_saf_id}_{bam_list_id}.tsv.summary"
    log:
              "out/{tool}{extra}_{gtf_saf_id}_{bam_list_id}.log"
    benchmark:
              "out/{tool}{extra}_{gtf_saf_id}_{bam_list_id}.benchmark.tsv"
    params:
        extra = params_extra
    conda:
        "../envs/subread.yaml"
    wildcard_constraints:
        tool="htseq/count",
        bam_list="bam-[a-zA-Z0-9-]+",
        gtf_saf_id="[a-zA-Z0-9-]+"
    threads:
        MAX_THREADS
    shell:
        """
        htseq-count -n {threads} {params.extra} {input.bam} {input.gtf} > {output.tsv} 2> {log}
        """

