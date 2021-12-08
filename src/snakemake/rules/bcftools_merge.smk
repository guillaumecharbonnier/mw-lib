rule bcftools_merge:
    """
    Test:
    """
    input:
        vcf_list  = lambda wildcards: eval(mwconf['ids'][wildcards.vcf_list_id]),
    output:
        "out/{tool}{extra}/{vcf_list_id}.vcf"
    params:
        extra = params_extra
    wildcard_constraints:
        extin="vcf|bcf",
        tool="bcftools/merge"
    conda:
        "../../../../mw-lib/src/snakemake/envs/bcftools.yaml"
    shell:
        """
        bcftools merge {params.extra} {input} > {output}
        """

