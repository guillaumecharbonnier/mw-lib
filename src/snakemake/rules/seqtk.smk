rule seqtk_seq_fastq_to_fasta:
    """
    Created:
        2017-04-12 17:26:22
    Aim:
        Seems to be the fastest way to convert fastq to fasta.
        https://www.biostars.org/p/85929/
    Doc:
        https://github.com/lh3/seqtk
    Test:
        out/seqtk/seq_A/sra-tools/fastq-dump_se/SRR1202037.fasta
    """
    input:
        fastq="out/{filler}.fastq"
    output:
        fasta="out/seqtk/seq_A/{filler}.fasta"
    conda:
        "../envs/seqtk.yaml"
    shell:
        "seqtk seq -A {input.fastq} > {output.fasta}"

rule seqtk_sample:
    """
    Aim:
        Seems to be the fastest way to convert fastq to fasta.
        https://www.biostars.org/p/85929/
    Doc:
        https://github.com/lh3/seqtk
    Test:
        out/seqtk/sample_100000/sra-tools/fastq-dump_se/SRR1202037.fastq.gz
    """
    input:
        "out/{filler}.{ext}"
    output:
        "out/{tool}{extra}_{n_or_frac}/{filler}.{ext}"
    params:
        extra = params_extra,
        gz = lambda wildcards: "| gzip -c" if wildcards.ext == 'fastq.gz' else ""
    wildcard_constraints:
        tool = "seqtk/sample",
        ext = "fastq|fastq.gz|fasta|fasta.gz",
        n_or_frac = "[0-9]+|0.[0-9]+"
    conda:
        "../envs/seqtk.yaml"
    shell:
        "seqtk sample {params.extra} {input} {wildcards.n_or_frac} {params.gz} > {output}"

