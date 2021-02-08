rule deepvariant_fa_bam_dev_only:
    """
    Created:
        2020-01-22 22:18:19
    Warning:
        Deepvariant is not deployed easily without Docker so this is just a rule to test it locally.
    Doc:
        https://github.com/google/deepvariant/blob/r0.9/docs/deepvariant-quick-start.md
    Note:
        Wrapper script is not currently available in bioconda package so git clone deepvariant git repository first.
    Test:
        out/deepvariant/_fa-genome-hg19-main-chr/samtools/index/abra2/_--single_fa-genome-hg19-main-chr/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_se_fa-genome-hg19-main-chr/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_H3K27ac.vcf.gz
    """
    input:
        script="../deepvariant/scripts/run_deepvariant.py",
        bam="out/{filler}.bam",
        fa= lambda wildcards: eval(mwconf['ids'][wildcards.fa_genome_id])
    output:
        vcf="out/{tool}{extra}_{fa_genome_id}/{filler}.vcf.gz",
        gvcf="out/{tool}{extra}_{fa_genome_id}/{filler}.g.vcf.gz"
    log:
            "out/{tool}{extra}_{fa_genome_id}/{filler}.log"
    benchmark:
            "out/{tool}{extra}_{fa_genome_id}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    #conda:
    #    "../envs/deepvariant.yaml"
    wildcard_constraints:
        tool="deepvariant/"
    threads:
        1
    shell:
        """
        mkdir -p `dirname {output.vcf}`
        BIN_VERSION="0.9.0"
        sudo docker run -v "`pwd`":"/mw/" \
          google/deepvariant:"${{BIN_VERSION}}" \
          /opt/deepvariant/bin/run_deepvariant \
          --model_type=WGS \
          --ref=/mw/{input.fa} \
          --reads=/mw/{input.bam} \
          --regions "chr1:4,7670,000-4,780,000" \
          --output_vcf=/mw/{output.vcf} \
          --output_gvcf=/mw/{output.gvcf} \
          --num_shards={threads}
        """

