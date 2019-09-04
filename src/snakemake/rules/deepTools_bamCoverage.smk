rule deepTools_bamCoverage_extra:
    """
    Modified:
        2016-10-17 10h27 - This rule should include the cases of two commented following rules, making them obsolete.
        2017-09-27 17:50:41 - Remove set values for argument --minMappingQuality, --binSize and --normalizeUsingRPKM to follow new paradigm.
        2018-12-06 10:57:35 - Changed to params_extra paradigm.
    Aim:
        This rule use DeepTools bamCoverage to transfrom bam to bigWig.
    Test:
        out/deepTools/bamCoverage/bowtie2/filler/sra-tools/fastq-dump/SRR1202037.bw
        out/deepTools/bamCoverage_--filterRNAstrand_reverse_--binSize_20_--minMappingQuality_0_--normalizeUsing_RPKM/star/se_hg19/sickle/se_-t_sanger_-q_20/gunzip/merge_lanes_nextseq500_se_raw/inp/fastq/run145/S001375_K562_IFN_R2-18036/K562-IFN-R2_S7.bw
    """
    input:
        bam="out/{filler}.bam",
        bai="out/{filler}.bam.bai"
    output:
        bw="out/{tool}{extra}/{filler}.bw"
    log:
        "out/{tool}{extra}/{filler}.log"
    benchmark:
        "out/{tool}{extra}/{filler}.benchmark.tsv"
    wildcard_constraints:
        tool="deepTools/bamCoverage"
    params:
        extra = params_extra
    conda:
        "../envs/deeptools.yaml"
    threads:
        16
    shell:
        "bamCoverage --bam {input.bam} --numberOfProcessors {threads} {params.extra} -o {output.bw} &> {log}"

#####################
# ONLY LEGACY BELOW #
#####################

rule deepTools_bamCoverage_binSize_minMappingQuality_legacy:
    """
    Modified:
        2017-01-24 17h59
        2017-03-08 10:22:41 - Renamed rule to new paradigm.
        2017-09-27 17:53:05 - Updated to new paradigm with '-' between argument name and value. Also removed --normalizeUsingRPKM
    Aim:
        This rule use DeepTools bamCoverage to transfrom bam to bigWig.
    Test:
    """
    input:
        bamCoverage="opt/miniconda/envs/deeptools/bin/bamCoverage",
        bam="out/{filler}.bam",
        bai="out/{filler}.bam.bai"
    output:
        bw="out/deepTools/bamCoverage_binSize-{binSize}_minMappingQuality-{minMappingQuality}/{filler}.bw"
    #conda:
    #    "../envs/deeptools.yaml"
    wildcard_constraints:
        binSize="[0-9]+",
        minMappingQuality="[0-9]+"
    threads:
        16
    shell:
        """
        {input.bamCoverage} \
            --bam {input.bam} \
            --numberOfProcessors {threads} \
            --binSize {wildcards.binSize} \
            --minMappingQuality {wildcards.minMappingQuality} \
            -o {output.bw}
        """

rule deepTools_bamCoverage_binSize_minMappingQuality_normalizeUsingRPKM_legacy:
    """
    Created:
        2017-04-10 15:18:14
    Aim:
        This rule use DeepTools bamCoverage to transfrom bam to bigWig.
    """
    input:
        bamCoverage="opt/miniconda/envs/deeptools/bin/bamCoverage",
        bam="out/{filler}.bam",
        bai="out/{filler}.bam.bai"
    output:
        bw="out/deepTools/bamCoverage_binSize-{binSize}_minMappingQuality-{minMappingQuality}_normalizeUsingRPKM/{filler}.bw"
    conda:
        "../envs/deeptools.yaml"
    wildcard_constraints:
        binSize="[0-9]+",
        minMappingQuality="[0-9]+"
    threads:
        16
    shell:
        """
        {input.bamCoverage} \
            --bam {input.bam} \
            --numberOfProcessors {threads} \
            --binSize {wildcards.binSize} \
            --minMappingQuality {wildcards.minMappingQuality} \
            --normalizeUsingRPKM \
            --outFileName {output.bw}
        """

rule deepTools_bamCoverage_binSize_extendReads_legacy:
    """
    Created:
        2018-03-05 15:39:13
    Aim:
        This rule use DeepTools bamCoverage to transfrom bam to bigWig.
    Note:
    Test:
    """
    input:
        bamCoverage="opt/miniconda/envs/deeptools/bin/bamCoverage",
        bam="out/{filler}.bam",
        bai="out/{filler}.bam.bai"
    output:
        bw="out/deepTools/bamCoverage_binSize-{binSize}_extendReads-{extendReads}/{filler}.bw"
    wildcard_constraints:
        binSize="[0-9]+",
        extendReads="|[0-9]+", # empty argument is for paired-end data.
    threads:
        16
    shell:
        """
        {input.bamCoverage}\
            --bam {input.bam}\
            --numberOfProcessors {threads}\
            --binSize {wildcards.binSize}\
            --extendReads {wildcards.extendReads}\
            --outFileName {output.bw}
        """

rule deepTools_bamCoverage_binSize_normalizeUsing_extendReads_legacy:
    """
    Created:
        2017-10-30 09:27:14
    Modified:
        2018-06-27 18:49:57 - --normalizeUsingRPKM does not work anymore.
    Aim:
        This rule use DeepTools bamCoverage to transfrom bam to bigWig.
    Note:
    Test:
    """
    input:
        bamCoverage="opt/miniconda/envs/deeptools/bin/bamCoverage",
        bam="out/{filler}.bam",
        bai="out/{filler}.bam.bai"
    output:
        bw="out/deepTools/bamCoverage_binSize-{binSize}_normalizeUsing-{normalizeUsing}_extendReads-{extendReads}/{filler}.bw"
    wildcard_constraints:
        binSize="[0-9]+",
        normalizeUsing="RPKM|CPM|BPM|RPGC|None",
        extendReads="|[0-9]+", # empty argument is for paired-end data.
    threads:
        16
    shell:
        """
        {input.bamCoverage} \
            --bam {input.bam} \
            --numberOfProcessors {threads} \
            --binSize {wildcards.binSize} \
            --normalizeUsing {wildcards.normalizeUsing} \
            --extendReads {wildcards.extendReads} \
            --outFileName {output.bw}
        """

rule deepTools_bamCoverage_binSize_minMappingQuality_normalizeUsingRPKM_extendReads_legacy:
    """
    Created:
        2017-05-08 14:59:15 - extendReads is useful for paired-end data.
    Aim:
        This rule use DeepTools bamCoverage to transfrom bam to bigWig.
    Note:
        extendReads with empty value seems bugged in some deepTools versions.
    Test:
        out/deepTools/bamCoverage_binSize-10_minMappingQuality-0_normalizeUsingRPKM_extendReads-/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run113_run125_run126/PSK-SC-WT.bw
    """
    input:
        bamCoverage="opt/miniconda/envs/deeptools/bin/bamCoverage",
        bam="out/{filler}.bam",
        bai="out/{filler}.bam.bai"
    output:
        bw="out/deepTools/bamCoverage_binSize-{binSize}_minMappingQuality-{minMappingQuality}_normalizeUsingRPKM_extendReads-{extendReads}/{filler}.bw"
    #conda:
    #    "../envs/deeptools.yaml"
    wildcard_constraints:
        binSize="[0-9]+",
        minMappingQuality="[0-9]+",
        extendReads="|[0-9]+", # empty argument is for paired-end data.
    threads:
        16
    shell:
        """
        {input.bamCoverage} \
            --bam {input.bam} \
            --numberOfProcessors {threads} \
            --binSize {wildcards.binSize} \
            --minMappingQuality {wildcards.minMappingQuality} \
            --normalizeUsingRPKM \
            --extendReads {wildcards.extendReads} \
            --outFileName {output.bw}
        """

rule deepTools_version_bamCoverage_binSize_minMappingQuality_normalizeUsingRPKM_extendReads_legacy:
    """
    Created:
        2017-05-08 14:59:15 - extendReads is useful for paired-end data.
    Aim:
        This rule use DeepTools bamCoverage to transfrom bam to bigWig.
    Note:
        extendReads with empty value seems bugged in some deepTools versions.
    Test:
        version 2.5.1 does not work since reinstallation on 2017-09-26 14:36:47. out/deepTools/2_5_1_bamCoverage_binSize-10_minMappingQuality-0_normalizeUsingRPKM_extendReads-/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run113_run125_run126/PSK-SC-WT.bw
        version 2.5.3 does.
    """
    input:
        bamCoverage="opt/miniconda/envs/deeptools_{version}/bin/bamCoverage",
        bam="out/{filler}.bam",
        bai="out/{filler}.bam.bai"
    output:
        bw="out/deepTools/{version}_bamCoverage_binSize-{binSize}_minMappingQuality-{minMappingQuality}_normalizeUsingRPKM_extendReads-{extendReads}/{filler}.bw"
    #conda:
    #    "../envs/deeptools.yaml"
    wildcard_constraints:
        version="[0-9]_[0-9]_[0-9]",
        binSize="[0-9]+",
        minMappingQuality="[0-9]+",
        extendReads="|[0-9]+", # empty argument is for paired-end data.
    threads:
        16
    shell:
        """
        export PATH="opt/miniconda/envs/deeptools_{wildcards.version}/bin:$PATH"

        {input.bamCoverage} \
            --bam {input.bam} \
            --numberOfProcessors {threads} \
            --binSize {wildcards.binSize} \
            --minMappingQuality {wildcards.minMappingQuality} \
            --normalizeUsingRPKM \
            --extendReads {wildcards.extendReads} \
            --outFileName {output.bw}
        """

rule deepTools_bamCoverage_filterRNAstrand_binSize_minMappingQuality_normalizeUsingRPKM_legacy:
    """
    Created:
        2017-03-28 10:39:21 - Forked to produce forward and reverse coverage for stranded RNA-Seq
    Aim:
        This rule use DeepTools bamCoverage to transfrom bam to bigWig.
    Test:
        out/deepTools/bamCoverage_filterRNAstrand-forward_binSize-20_minMappingQuality-0_normalizeUsingRPKM/star/pe_GRCm38_outFilterMultimapNmax-1000/sickle/pe_-t_sanger_-q_20/merge_lanes/run176/RNA-P-H2AL2-WT-Rep1.bw
    """
    input:
        bamCoverage="opt/miniconda/envs/deeptools/bin/bamCoverage",
        bam="out/{filler}.bam",
        bai="out/{filler}.bam.bai"
    output:
        bw="out/deepTools/bamCoverage_filterRNAstrand-{filterRNAstrand}_binSize-{binSize}_minMappingQuality-{minMappingQuality}_normalizeUsingRPKM/{filler}.bw"
    #conda:
    #    "../envs/deeptools.yaml"
    wildcard_constraints:
        filterRNAstrand="forward|reverse",
        binSize="[0-9]+",
        minMappingQuality="[0-9]+"
    threads:
        16
    shell:
        """
        {input.bamCoverage} \
            --bam {input.bam} \
            --numberOfProcessors {threads} \
            --filterRNAstrand {wildcards.filterRNAstrand} \
            --binSize {wildcards.binSize} \
            --minMappingQuality {wildcards.minMappingQuality} \
            --normalizeUsingRPKM \
            --outFileName {output.bw}
        """


#rule deepTools_bamCoverage_from_sam_to_bam:
#    """
#    Modified: 2016-09-28 -
#    This rule use DeepTools bamCoverage to transfrom bam to bigWig.
#    """
#    input:
#        bamCoverage="opt/miniconda/envs/deeptools/bin/bamCoverage",
#        bam="out/samtools/sam_to_bam/bowtie2/{runtype}/{index}/{exp}/{sample}.bam",
#        bai="out/samtools/sam_to_bam/bowtie2/{runtype}/{index}/{exp}/{sample}.bam.bai",
#    output:
#        bw="out/deepTools/bamCoverage/samtools/sam_to_bam/bowtie2/{runtype}/{index}/{exp}/{sample}.bw"
#    threads: 16
#    shell:"""
#    {input.bamCoverage} --bam {input.bam} \
#    --numberOfProcessors 16 --binSize 1 \
#    --minMappingQuality 30 \
#    --normalizeUsingRPKM -o {output.bw}
#    """
#
#rule deepTools_bamCoverage_from_star:
#    """
#    Modified: 2016-09-28 -
#    This rule use DeepTools bamCoverage to transfrom bam to bigWig.
#
#    out/star/se/mm9/GSE49622/SRR948773.bam
#    """
#    input:
#        bamCoverage="opt/miniconda/envs/deeptools/bin/bamCoverage",
#        bam="out/star/{runtype}/{index}/{exp}/{sample}.bam",
#        bai="out/star/{runtype}/{index}/{exp}/{sample}.bam.bai",
#    output:
#        bw="out/deepTools/bamCoverage/star/{runtype}/{index}/{exp}/{sample}.bw"
#    threads: 16
#    shell:"""
#    {input.bamCoverage} --bam {input.bam} \
#    --numberOfProcessors 16 --binSize 1 \
#    --minMappingQuality 30 \
#    --normalizeUsingRPKM -o {output.bw}
#    """

rule deepTools_bamCoverage_blacklisted_legacy:
    """
    Modified: 2016-04-25 Added minimum quality and removed extension of reads to 140pb.
    This rule use DeepTools bamCoverage to transfrom bam to bigWig.
    """
    input:
        bamCoverage="opt/miniconda/envs/deeptools/bin/bamCoverage",\
        bam="out/bam/{index}/{exp}/{sample}.bam",\
        bed="annotation/input/feature/{index}/blacklist.bed"
    output:
        bw="out/bw/bamCoverage_blacklisted/{index}/{exp}/{sample}.bw"
    threads: 16
    shell:"""
    {input.bamCoverage} --bam {input.bam} \
        --numberOfProcessors 16 --binSize 1 \
        --minMappingQuality 30 \
        --blackListFileName {input.bed} \
        --normalizeUsingRPKM -o {output.bw}
    """

rule deepTools_bamCoverage_MNase_mode_legacy:
    """
    Created: 2016-04-22 23h43 - Imported from "deepTools_bamCoverage"
    """
    input:
        bamCoverage="opt/miniconda/envs/deeptools/bin/bamCoverage",\
        bam="out/bam/{index}/{exp}/{sample}.bam"
    output: bw="out/bw/bamCoverage_MNase_mode/{index}/{exp}/{sample}.bw"
    threads: 16
    shell:"""
        {input.bamCoverage} --bam {input.bam} \
        --numberOfProcessors 16 --binSize 1 \
        --normalizeUsingRPKM -o {output.bw} \
        --MNase
    """

rule deepTools_bamCoverage_custom_MNase_mode:
    """
    Created: 2016-04-22 23h43 - Imported from "deepTools_bamCoverage"
    """
    input:
        bamCoverage="opt/miniconda/envs/deeptools/bin/bamCoverage",
        bam="out/bam/{index}/{exp}/{sample}.bam"
    output:
        bw="out/bw/bamCoverage_custom_MNase_mode/{index}/{exp}/{sample}.bw"
    threads: 16
    shell:"""
    {input.bamCoverage} --bam {input.bam} \
        --numberOfProcessors 16 --binSize 1 \
        --normalizeUsingRPKM -o {output.bw} \
        --centerReads --extendReads 3
    """

rule deepTools_bamCoverage_legacy:
    """
    Modified:
        2016-03-15 13h14 - Updated bamCoverage.
    Modified:
        2016-04-25 - Added minimum quality and removed extension of reads to 140pb.
    Modified:
        2017-03-08 10:21:02 - Changed to legacy because of deprecated file patterns.
    Aim:
        This rule use DeepTools bamCoverage to transfrom bam to bigWig.
    """
    input:
        #bamCoverage="opt/miniconda/envs/deeptools/bin/bamCoverage",
        bam="out/bam/{index}/{exp}/{sample}.bam"
    output:
        bw="out/bw/bamCoverage/{index}/{exp}/{sample}.bw"
    conda:
        "../envs/deeptools.yaml"
    threads: 16
    shell:"""
    #{input.bamCoverage} --bam {input.bam} \
    bamCoverage --bam {input.bam} \
    --numberOfProcessors 16 --binSize 1 \
    --minMappingQuality 30 \
    --normalizeUsingRPKM -o {output.bw}
    """
