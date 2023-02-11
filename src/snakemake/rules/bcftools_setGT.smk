rule bcftools_setGT_somatic_to_missing:
    """
    bcftools +setGT test.vcf -- -t q -n . -i'FORMAT/GP>=0.99'
    
    Test:
        out/bcftools/setGT_somatic_to_missing/bcftools/fill-tags_-t_VAF/inp/inp/vcf/3186patients_complet.sorted.vcf.gz
    """
    input:
        "out/{filler}.{ext}"
    output:
        "out/{tool}/{filler}.{ext}"
    params:
        #extra = params_extra,
        O = lambda wildcards: "-Ob" if wildcards.ext == 'bcf' else "-Oz" if wildcards.ext == 'vcf.gz' else ""
    wildcard_constraints:
        extin="vcf|bcf|vcf.gz",
        tool="bcftools/setGT_somatic_to_missing"
    conda:
        "../../../../mw-lib/src/snakemake/envs/bcftools.yaml"
    shell:
        """
        bcftools +setGT {input} {params.O} -o {output} -- -t q -n . -i'FORMAT/VAF<=0.45&FORMAT/VAF>0'
        """
