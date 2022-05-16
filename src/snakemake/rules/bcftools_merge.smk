rule bcftools_merge:
    """
    Test:
        out/bcftools/merge/vcf-GRCh38-TALL-H3K27ac-H3K4me3.vcf
    """
    input:
        vcf_list  = lambda wildcards: eval(mwconf['ids'][wildcards.vcf_list_id]),
    output:
        "out/{tool}{extra}/{vcf_list_id}.vcf"
    log:
        "out/{tool}{extra}/{vcf_list_id}.log"
    benchmark:
        "out/{tool}{extra}/{vcf_list_id}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        extin="vcf|bcf",
        tool="bcftools/merge"
    conda:
        "../../../../mw-lib/src/snakemake/envs/bcftools.yaml"
    shell:
        """
        bcftools merge {params.extra} {input} > {output} 2> {log}
        """

