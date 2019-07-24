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
