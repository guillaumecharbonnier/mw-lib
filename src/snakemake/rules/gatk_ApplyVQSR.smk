rule gatk4_ApplyVQSR:
    """
    Created:
        2022-03-29 17:47:08
    Aim:
        Recalibrate SNP
    Doc:
        https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
        https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.4.0/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php
    Note:
    Test:
        out/gatk4/ApplyVQSR/pbgzip_tabix/sed/rm-chr-in-vcf/gunzip/to-stdout/ln/updir/mw-tall-pediac/inp/vcf/doublons_captures_diff.vcf.gz
        https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.22.vcf.bgz' --out 'gnomad.exomes.r2.1.1.sites.22.vcf.bgz
    """
    input:
        vcf = "out/{filler}.vcf.gz",
        tbi = "out/{filler}.vcf.gz.tbi",
        recal    = "out/gatk4/VariantRecalibrator/{filler}.recal",
        tranches = "out/gatk4/VariantRecalibrator/{filler}.tranches"
    output:
        vcf="out/{tool}{extra}/{filler}.vcf.gz"
    log:
        "out/{tool}{extra}/{filler}.log"
    benchmark:
        "out/{tool}{extra}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "gatk4/ApplyVQSR"
    threads:
        MAX_THREADS
    conda:
        "../envs/gatk4_vqsr.yaml"
    shell:
        """
         gatk ApplyVQSR \
            -V {input.vcf} \
            -O {output.vcf} \
            --truth-sensitivity-filter-level 99.0 \
            --tranches-file {input.tranches} \
            --recal-file {input.recal} \
            -mode SNP
        """


