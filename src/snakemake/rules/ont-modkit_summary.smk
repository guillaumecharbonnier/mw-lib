rule ont_modkit_summary_extra:
    """
    Created:
        2024-11-06 17:14:04
    Doc:
        https://nanoporetech.github.io/modkit/
    Note:
    Test:
        out/ont-modkit/summary_--no-sampling/ln/updir/mw-tall-rrbs/inp/nanopore/DOUYERE_diag_hs1.bed
    """
    input:
        bam="out/{filler}.bam",
        bai="out/{filler}.bam.bai"
    output:
        txt="out/{tool}{extra}_{fa_genome_id}/{filler}.txt"
    log:
            "out/{tool}{extra}_{fa_genome_id}/{filler}.log"
    benchmark:
            "out/{tool}{extra}_{fa_genome_id}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    conda:
        "../envs/ont-modkit.yaml"
    wildcard_constraints:
        tool="ont-modkit/summary"
    threads:
        MAX_THREADS
    shell:
        """
        modkit summary {params.extra} --threads {threads} --log {log}  {input.bam} {output.txt}
        """
