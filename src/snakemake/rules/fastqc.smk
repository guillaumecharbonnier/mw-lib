rule fastqc:
    """
    Created:
        2016-10-09-27 11:28:58
    Aim:
        Check read quality
    Doc:
        https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
    Test:
        out/fastqc/fastq/sra-tools/fastq-dump_se/SRR3126243_fastqc.html
    """
    input:
        "out/{filler}.{ext}"
    output:
        zip  = "out/{tool}{extra}{ext}/{filler}_fastqc.zip",
        html = "out/{tool}{extra}{ext}/{filler}_fastqc.html"
    log:
               "out/{tool}{extra}{ext}/{filler}.log"
    benchmark:
               "out/{tool}{extra}{ext}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    conda:
        "../envs/fastqc.yaml"
    wildcard_constraints:
        tool="fastqc/",
        ext="fastq|fastq.gz|bam|sam"
    threads:
        MAX_THREADS
    shell:
        "fastqc --threads {threads} --outdir `dirname {output.html}` {params.extra} {input} &> {log}"


