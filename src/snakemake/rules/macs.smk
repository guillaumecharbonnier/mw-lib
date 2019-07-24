"""
Warning:
    Rules below are just copy-paste of macs2 ones. Not adjusted to macs1 behaviour.
"""

rule macs_withctrl_extra:
    """
    Created:
        2019-02-18 15:27:59
    Doc:
        https://github.com/taoliu/MACS/blob/macs_v1/README.rst
    Aim:
        Call peaks with control/input
    Test:
        out/macs/withctrl/samtools/index/samtools/sort/samtools/view_sam_to_bam_-q_30/bowtie2/se_mm10/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR3126243_over_SRR3126242_peaks.bed
    """
    input:
        bam_chip  = "out/{filler}/{chip}.bam",
        bai_chip  = "out/{filler}/{chip}.bam.bai",
        bam_input = "out/{filler}/{input}.bam",
        bai_input = "out/{filler}/{input}.bam.bai"
    output:
        bed    = "out/{tool}{extra}/{filler}/{chip}_over_{input}_peaks.bed",
        xls    = "out/{tool}{extra}/{filler}/{chip}_over_{input}_peaks.xls"
    params:
        outdir = "out/{tool}{extra}/{filler}",
        extra  = params_extra
    wildcard_constraints:
        tool = "macs/withctrl"
    conda:
        "../envs/macs.yaml"
    wildcard_constraints:
        chip  = "[a-zA-Z0-9-_]+",
        input = "[a-zA-Z0-9-_]+"
    shell:
        """
        macs {params.extra}\
            --treatment {input.bam_chip}\
            --control {input.bam_input}\
            --name {wildcards.chip}_over_{wildcards.input}\
            --outdir {params.outdir}
        # Renaming narrowPeak or broadPeak to have only one output bed name for all variations of settings
        #TO_RENAME=`find {params.outdir} -name '*narrowPeak' -o -name '{wildcards.chip}_over_{wildcards.input}_peaks.*broadPeak'`
        #ln $TO_RENAME {output.bed}
        """

rule macs_noctrl_extra:
    """
    Created:
        2018-10-22 12:33:01
    Aim:
        Call peaks without control/input.
    Test:
        out/macs/noctrl/samtools/index/samtools/sort/samtools/view_sam_to_bam_-q_30/bowtie2/se_mm10/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR3126243_peaks.bed
    """
    input:
        bam_chip  = "out/{filler}/{chip}.bam",
        bai_chip  = "out/{filler}/{chip}.bam.bai",
    output:
        bed    = "out/{tool}{extra}/{filler}/{chip}_peaks.bed",
        xls    = "out/{tool}{extra}/{filler}/{chip}_peaks.xls"
    params:
        outdir = "out/{tool}{extra}/{filler}",
        extra  = params_extra
    wildcard_constraints:
        tool = "macs/noctrl"
    conda:
        "../envs/macs.yaml"
    wildcard_constraints:
        chip="[a-zA-Z0-9-_]+"
    shell:
        """
        cd {params.outdir}
        macs {params.extra}\
            --treatment {WDIR}/{input.bam_chip}\
            --name {wildcards.chip}

        #    --outdir {params.outdir}
        # Renaming narrowPeak or broadPeak to have only one output bed name for all variations of settings
        #TO_RENAME=`find {params.outdir} -name '*narrowPeak' -o -name '*broadPeak'`
        #ln $TO_RENAME {output.bed}
        """


