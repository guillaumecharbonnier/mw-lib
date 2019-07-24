rule cutadapt_single_end_extra:
    """
    Created:
        2019-02-01 18:23:30
    Aim:

    Note:
        --nextseq-trim=3'CUTOFF
        NextSeq-specific quality trimming (each read). Trims
        also dark cycles appearing as high-quality G bases
        (EXPERIMENTAL).
    Test:
        out/cutadapt/se_-q_30/sra-tools/fastq-dump_se/SRR1202037.fastq
    """
    input:
        "out/{filler}.{ext}"
    output:
        "out/{tool}{extra}/{filler}.{ext}"
    log:
        "out/{tool}{extra}/{filler}.{ext}.log"
    benchmark:
        "out/{tool}{extra}/{filler}.{ext}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="cutadapt/se",
        ext="(fasta|fa|fna|fastq|fq)(\.gz|\.bz2|\.xz|)"
    conda:
        "../envs/cutadapt.yaml"
    threads:
        16
    shell:
        "cutadapt {params.extra} -j {threads} -o {output} {input} &> {log}"

