rule ont_modkit_pileup_extra:
    """
    Created:
        2024-11-06 17:14:04
    Doc:
        https://nanoporetech.github.io/modkit/
    Note:
    Test:
        #out/ont-modkit/pileup_--cpg_fa-genome-hs1/ln/updir/mw-tall-rrbs/inp/nanopore/DOUYERE_diag_hs1.bedmethyl
        out/ont-modkit/pileup_--cpg_fa-genome-hs1/ln/updir/mw-tall-rrbs/inp/nanopore/DOUYERE_diag_hs1.bed
    """
    input:
        bam="out/{filler}.bam",
        bai="out/{filler}.bam.bai",
        fa= lambda wildcards: eval(mwconf['ids'][wildcards.fa_genome_id])
    output:
        bed="out/{tool}{extra}_{fa_genome_id}/{filler}.bed"
    log:
            "out/{tool}{extra}_{fa_genome_id}/{filler}.log"
    benchmark:
            "out/{tool}{extra}_{fa_genome_id}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    conda:
        "../envs/ont-modkit.yaml"
    wildcard_constraints:
        tool="ont-modkit/pileup"
    threads:
        MAX_THREADS
    shell:
        """
        modkit pileup {params.extra} --threads {threads} --ref {input.fa} --log {log} {input.bam} {output.bed}
        """
