rule samtools_index:
    """
    Modified:
        2017-04-11 14:52:20 - Exported to its own file.
        2018-01-09 13:35:50 - Adding the 'out' prefix to prevent ambiguity issues with rule ln_root_dir_to_out. NOT TESTED. See ruleorder debug in includerules.
        2019-01-28 20:01:25 - Adding the 'samtools/index' prefix to prevent ambiguity issues with crossmap_bam rule. Should break most of the old code but great improvement to behaviour. Samtools index is now a rule like any other in MetaWorkflow.
    Aim:
        Index bam file for use with tools like genome browsers.
    Test:

    """
    input:
        bam="out/{filler}.bam",
    output:
        bam="out/samtools/index/{filler}.bam",
        bai="out/samtools/index/{filler}.bam.bai"
    log:
            "out/samtools/index/{filler}.log"
    benchmark:
            "out/samtools/index/{filler}.benchmark.tsv"
    conda:
        "../envs/samtools.yaml"
    envmodules:
        "samtools/1.14"
    threads:
        1
    shell:
        """
        (ln -srf {input.bam} {output.bam}
        samtools index {output.bam}) &> {log}
        """

ruleorder: crossmap_bam > samtools_index_legacy
ruleorder: ln_srf_parent_dir > samtools_index_legacy
ruleorder: ln_alias > samtools_index_legacy
ruleorder: samtools_index > samtools_index_legacy
ruleorder: samtools_sam_to_bam_bai_extra > samtools_index_legacy
rule samtools_index_legacy:
    """
    Modified:
        2017-04-11 14:52:20 - Exported to its own file.
        2018-01-09 13:35:50 - Adding the 'out' prefix to prevent ambiguity issues with rule ln_root_dir_to_out. NOT TESTED. See ruleorder debug in includerules.
        2019-01-28 20:01:25 - Adding the 'samtools/index' prefix to prevent ambiguity issues with crossmap_bam rule. Should break most of the old code but great improvement to behaviour. Samtools index is now a rule like any other in MetaWorkflow.
    Aim:
        Index bam file for use with tools like genome browsers.
    Test:
    Todo:
        Replace this rule everywhere with the new samtools_index one.
    """
    input:
        bam="out/{filler}.bam",
    output:
        bai="out/{filler}.bam.bai"
    log:
        "out/{filler}.bam.bai.log"
    benchmark:
        "out/{filler}.bam.bai.benchmark.tsv"
    conda:
        "../envs/samtools.yaml"
    threads:
        1
    shell:
        "samtools index {input.bam} &> {log}"
