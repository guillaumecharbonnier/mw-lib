rule gatk4_VariantRecalibrator:
    """
    Created:
        2022-03-29 17:47:08
    Aim:
        Recalibrate SNP
    Doc:
        https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
    Test:
        out/gatk4/VariantRecalibrator/pbgzip_tabix//sed/rm-chr-in-vcf/gunzip/to-stdout/ln/updir/mw-tall-pediac/inp/vcf/doublons_captures_diff.recal
    """
    input:
        vcf="out/{filler}.vcf.gz",
        tbi="out/{filler}.vcf.gz.tbi",
        # hapmap = "out/wget/https/storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz",
        # omni   = "out/wget/https/storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz",
        # phase1 = "out/wget/https/storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        # dbsnp  = "out/wget/https/storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf"
        hapmap = "out/pbgzip_tabix/gunzip/to-stdout/wget/https/storage.googleapis.com/gatk-legacy-bundles/b37/hapmap_3.3.b37.vcf.gz",
        omni   = "out/pbgzip_tabix/gunzip/to-stdout/wget/https/storage.googleapis.com/gatk-legacy-bundles/b37/1000G_omni2.5.b37.vcf.gz",
        phase1 = "out/pbgzip_tabix/gunzip/to-stdout/wget/https/storage.googleapis.com/gatk-legacy-bundles/b37/1000G_phase1.indels.b37.vcf.gz",
        dbsnp  = "out/pbgzip_tabix/gunzip/to-stdout/wget/https/storage.googleapis.com/gatk-legacy-bundles/b37/dbsnp_138.b37.vcf.gz"
    output:
        recal    = "out/{tool}{extra}/{filler}.recal",
        tranches = "out/{tool}{extra}/{filler}.tranches",
        rscript  = "out/{tool}{extra}/{filler}.R",
    log:
        "out/{tool}{extra}/{filler}.log"
    benchmark:
        "out/{tool}{extra}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "gatk4/VariantRecalibrator"
    threads:
        MAX_THREADS
    conda:
        "../envs/gatk4_vqsr.yaml"
    shell:
        """
        gatk --java-options "-Xmx3g -Xms3g" VariantRecalibrator \
            -V {input.vcf} \
            --trust-all-polymorphic \
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
            -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
            -mode SNP \
            --max-gaussians 6 \
            -resource:hapmap,known=false,training=true,truth=true,prior=15 {input.hapmap} \
            -resource:omni,known=false,training=true,truth=true,prior=12 {input.omni} \
            -resource:1000G,known=false,training=true,truth=false,prior=10 {input.phase1} \
            -resource:dbsnp,known=true,training=false,truth=false,prior=7 {input.dbsnp} \
            -O {output.recal} \
            --tranches-file {output.tranches} \
            --rscript-file {output.rscript}
        """


