rule bcftools_filltags:
    """
    Test:
        out/bcftools/fill-tags_-t_VAF/inp/inp/vcf/3186patients_complet.sorted.vcf.gz
    """
    input:
        "out/{filler}.{ext}"
    output:
        "out/{tool}_{extra}/{filler}.{ext}"
    params:
        extra = params_extra,
        O = lambda wildcards: "-Ob" if wildcards.ext == 'bcf' else "-Oz" if wildcards.ext == 'vcf.gz' else ""
    wildcard_constraints:
        extin="vcf|bcf|vcf.gz",
        tool="bcftools/fill-tags"
    conda:
        "../../../../mw-lib/src/snakemake/envs/bcftools.yaml"
    shell:
        """
        bcftools +fill-tags {input} {params.O} -o {output} -- {params.extra}
        """
