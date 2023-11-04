rule bcftools_setGT_tests:
    input:
        "out/bcftools/setGT_low-qual-to-missing/bcftools/fill-tags_-t_VAF/ln/updir/mw-tall-pediac-data/inp/vcf/3186patients_complet.sorted.vcf.gz",
        "out/bcftools/setGT_somatic-to-missing/bcftools/fill-tags_-t_VAF/ln/updir/mw-tall-pediac-data/inp/vcf/3186patients_complet.sorted.vcf.gz",
        "out/bcftools/setGT_low-qual-to-missing/bcftools/setGT_somatic-to-missing/bcftools/fill-tags_-t_VAF/ln/updir/mw-tall-pediac-data/inp/vcf/3186patients_complet.sorted.vcf.gz",
        "out/bcftools/setGT_low-qual-or-somatic-to-missing/bcftools/fill-tags_-t_VAF/ln/updir/mw-tall-pediac-data/inp/vcf/3186patients_complet.sorted.vcf.gz"
    conda:
        "../../../../mw-lib/src/snakemake/envs/bcftools.yaml"
    shell:
        """
        for VCF in {input}; do
            echo $VCF
            bcftools view -H $VCF | md5sum
        done
        """

rule bcftools_setGT_extra:
    """
    Aim:
        Edit genotypes in a VCF file.
    Note:
        Refer to http://samtools.github.io/bcftools/howtos/filtering.html
        for more information on the -i option.
    """
    input:
        "out/{filler}.{ext}"
    output:
        "out/{tool}{extra}/{filler}.{ext}"
    params:
        extra = params_extra,
        O = lambda wildcards: "-Ob" if wildcards.ext == 'bcf' else "-Oz" if wildcards.ext == 'vcf.gz' else ""
    wildcard_constraints:
        extin = "vcf|bcf|vcf.gz",
        tool = "bcftools/setGT"
    conda:
        "../../../../mw-lib/src/snakemake/envs/bcftools.yaml"
    shell:
        """
        bcftools +setGT {input} {params.O} -o {output} -- {params.extra}
        """


# rule bcftools_setGT_somatic_to_missing:
#     """
#     bcftools +setGT test.vcf -- -t q -n . -i'FORMAT/GP>=0.99'
#     """
#     input:
#         "out/{filler}.{ext}"
#     output:
#         "out/{tool}/{filler}.{ext}"
#     params:
#         #extra = params_extra,
#         O = lambda wildcards: "-Ob" if wildcards.ext == 'bcf' else "-Oz" if wildcards.ext == 'vcf.gz' else ""
#     wildcard_constraints:
#         extin="vcf|bcf|vcf.gz",
#         tool="bcftools/setGT_somatic_to_missing"
#     conda:
#         "../../../../mw-lib/src/snakemake/envs/bcftools.yaml"
#     shell:
#         """
#         bcftools +setGT {input} {params.O} -o {output} -- -t q -n . -i'FMT/VAF<=0.45&FMT/VAF>0'
#         """

# rule bcftools_setGT_low_qual_to_missing:
#     """
#     "GQ<20"
#     "DP<10"
#     Note: 
#         Refer to http://samtools.github.io/bcftools/howtos/filtering.html
#         for more information on the -i option.
#     Test:
#         out/bcftools/setGT_low_qual_to_missing/bcftools/fill-tags_-t_VAF/ln/updir/mw-tall-pediac-data/inp/vcf/3186patients_complet.sorted.vcf.gz
#     """
#     input:
#         "out/{filler}.{ext}"
#     output:
#         "out/{tool}/{filler}.{ext}"
#     params:
#         #extra = params_extra,
#         O = lambda wildcards: "-Ob" if wildcards.ext == 'bcf' else "-Oz" if wildcards.ext == 'vcf.gz' else ""
#     wildcard_constraints:
#         extin="vcf|bcf|vcf.gz",
#         tool="bcftools/setGT_low_qual_to_missing"
#     conda:
#         "../../../../mw-lib/src/snakemake/envs/bcftools.yaml"
#     shell:
#         """
#         bcftools +setGT {input} {params.O} -o {output} -- -t q -n . -i'FMT/GQ<20|FMT/DP<10'
#         """

# rule bcftools_setGT_low_qual_or_somatic_to_missing:
#     """
#     bcftools +setGT test.vcf -- -t q -n . -i'FORMAT/GP>=0.99'
#     Note: 
#         Refer to http://samtools.github.io/bcftools/howtos/filtering.html
#         for more information on the -i option.
#     Test:
#         out/bcftools/setGT_low_qual_or_somatic_to_missing/bcftools/fill-tags_-t_VAF/ln/updir/mw-tall-pediac-data/inp/vcf/3186patients_complet.sorted.vcf.gz
#     """
#     input:
#         "out/{filler}.{ext}"
#     output:
#         "out/{tool}/{filler}.{ext}"
#     params:
#         #extra = params_extra,
#         O = lambda wildcards: "-Ob" if wildcards.ext == 'bcf' else "-Oz" if wildcards.ext == 'vcf.gz' else ""
#     wildcard_constraints:
#         extin="vcf|bcf|vcf.gz",
#         tool="bcftools/setGT_low_qual_or_somatic_to_missing"
#     conda:
#         "../../../../mw-lib/src/snakemake/envs/bcftools.yaml"
#     shell:
#         """
#         bcftools +setGT {input} {params.O} -o {output} -- -t q -n . -i'FMT/VAF<=0.45&FMT/VAF>0|FMT/GQ<20|FMT/DP<10'
#         """
