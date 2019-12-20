rule breakdancer:
    """
    Doc:
        https://github.com/genome/breakdancer
    """
    input:
        bam="out/{filler}.bam",
        fa= lambda wildcards: eval(config['ids'][wildcards.fa_genome_id])
    output:
        vcf="out/{tool}{extra}_{fa_genome_id}/{filler}.vcf"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="breakdancer/max"
    conda:
        "../envs/breakdancer.yaml"
    shell:
        """
        platypus callVariants --bamFiles=input.bam --refFile=reference.fa --output=variant_calls.vcf
        """


