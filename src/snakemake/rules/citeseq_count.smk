rule cite_seq_count:
    """
    Add test lines
    Need to add cells number
    """
    input:
        fwd ="out/{filler}_1.fastq.gz",
        rev ="out/{filler}_2.fastq.gz",
        csv="out/{tool}/{extra}/{filler}/tags.csv"
    output:
        yaml="out/{tool}/{extra}/{filler}/run_report.yaml"
    conda:
        "../envs/cite_seq_count.yaml"
    params:
        extra = params_extra
    threads:
        8
    wildcard_constraints:
        tool="cite-seq_count"
    log:
        "out/{tool}/{extra}/{filler}/cite_seq_count.log"
    shell:'''
        (
        echo {input.fwd}
        echo {input.rev}
        echo {input.csv}
        echo {output.yaml}
        CITE-seq-Count {params.extra} -R1 {input.fwd} -R2 {input.rev} -t {input.csv} ) > {log}
        '''

