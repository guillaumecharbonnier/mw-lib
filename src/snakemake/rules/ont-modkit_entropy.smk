rule ont_modkit_entropy_extra:
    """
    Created:
        2024-11-06 17:14:04
    Doc:
        https://nanoporetech.github.io/modkit/
    Note:
    Test:
        out/ont-modkit/entropy_--cpg_fa-genome-hs1/ln/updir/mw-tall-rrbs/inp/nanopore/DOUYERE_diag_hs1.bed
        out/ont-modkit/entropy_--cpg_fa-genome-hs1/samtools/index/ln/updir/mw-tall-rrbs/inp/WGS/hs1_bams/ACHITE_diag_hs1.bed

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
    priority:
        5
    wildcard_constraints:
        tool="ont-modkit/entropy"
    threads:
        MAX_THREADS
    shell:
        """
        modkit entropy {params.extra} --threads {threads} --ref {input.fa} --log-filepath {log} -o {output.bed} --in-bam {input.bam} 
        """
