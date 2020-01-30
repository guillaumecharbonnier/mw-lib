rule stringtie:
    """
    Created:
        2020-01-30 00:56:54
    Doc:
        https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
    Note:
        Stopping rule writing because the tool appears not suitable for my current need and I would rather use RNASPades.
        RULE NOT COMPLETE
    Test:
    """
    input:
        bam="out/{filler}.bam",
        fa= lambda wildcards: eval(config['ids'][wildcards.fa_genome_id])
    output:
        vcf="out/{tool}{extra}_{fa_genome_id}/{filler}.vcf"
    log:
            "out/{tool}{extra}_{fa_genome_id}/{filler}.log"
    benchmark:
            "out/{tool}{extra}_{fa_genome_id}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    conda:
        "../envs/stringtie.yaml"
    wildcard_constraints:
        tool="stringtie/"
    threads:
        1
    shell:
        """
        #stringtie <aligned_reads.bam> [options]
        #freebayes {params.extra} -f {input.fa} {input.bam} 1> {output.vcf} 2> {log}
        """

