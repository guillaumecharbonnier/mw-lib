rule gatk4_ASEReadCounter:
    """
    Created:
        2022-03-29 17:47:08
    Aim:
        
    Note:
    Test:
        out/gatk4/ASEReadCounter_vcf-GRCh38-TP73-alleles/ln/alias/sst/all_samples/GRCh38/bam/1310_RNA_XT_HS2.tsv
    """
    input:
        bam = "out/{filler}.bam",
        vcf = lambda wildcards: eval(mwconf['ids'][wildcards.vcf_id])
    output:
        tsv="out/{tool}{extra}_{vcf_id}/{filler}.tsv"
    log:
        "out/{tool}{extra}_{vcf_id}/{filler}.log"
    benchmark:
        "out/{tool}{extra}_{vcf_id}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "gatk4/ASEReadCounter"
    threads:
        MAX_THREADS
    conda:
        # "../envs/gatk4_vqsr.yaml"
        "../envs/gatk4.yaml"
    shell:
        """
         gatk ASEReadCounter \
            --input {input.bam} \
            --variant {input.vcf} \
            -O {output.tsv}
        """


