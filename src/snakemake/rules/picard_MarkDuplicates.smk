ruleorder: picard_MarkDuplicates_RemoveDuplicates_AssumeSorted > picard_MarkDuplicates_extra

rule picard_MarkDuplicates_extra:
    """
    Created:
        2018-10-25 14:12:01
    Test:
        out/picard/MarkDuplicates_REMOVE_DUPLICATES=true/inp/bam/hg19_RNA-Seq_thymus/CD34plus.bam
    """
    input:
        bam="out/{filler}.bam",
        bai="out/{filler}.bam.bai"
    output:
        bam =     "out/{tool}{extra}/{filler}.bam",
        metrics = "out/{tool}{extra}/{filler}.metrics.txt"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "picard/MarkDuplicates"
    conda:
        "../envs/picard.yaml"
    shell:
        "picard MarkDuplicates I={input.bam} O={output.bam} M={output.metrics} {params.extra}"


#rule picard_MarkDuplicates:
#    """
#    Created:
#        2015
#    """
#    input:
#        picard="opt/miniconda/envs/picard/bin/picard",
#        samtools="opt/samtools-1.3.1/samtools",
#        bam="out/{id}.bam",
#        bai="out/{id}.bam.bai"
#    output:
#        bam =      "out/picard/MarkDuplicates/{id}.bam",
#        #bai =      "out/picard/MarkDuplicates/{id}.bam.bai",
#        metrics =  "out/picard/MarkDuplicates/{id}.metrics.txt",
#        flagstat = "out/picard/MarkDuplicates/{id}.flagstat.txt"
#    shell:
#        """
#        {input.picard} MarkDuplicates I={input.bam} O={output.bam} M={output.metrics} REMOVE_DUPLICATES=true
#        #{input.samtools} index {output.bam}
#        {input.samtools} flagstat {output.bam} > {output.flagstat}
#        """

rule picard_MarkDuplicates_RemoveDuplicates_AssumeSorted:
    """
    Created:
        2015

    java -Xmx2048m -jar picard/MarkDuplicates.jar INPUT=input.bam OUTPUT=output.bam METRICS_FILE=output.dup_metrics REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT

    Test:
       out/picard/MarkDuplicates_RemoveDuplicates-false_AssumeSorted-true/picard/SortSam_sortOrder-coordinate/samtools/view_bSh/bwa/samse_q-5_fa-GRCh38-Blueprint/gunzip/merge_lanes_nextseq500_se_raw/inp/fastq/run170/S001580_TH125_EC_K27ac-41046/TH125-EC-K27ac_S6.bam 
    """
    input:
        bam="out/{id}.bam"
    output:
        bam =     "out/picard/MarkDuplicates_RemoveDuplicates-{RemoveDuplicates}_AssumeSorted-{AssumeSorted}/{id}.bam",
        metrics = "out/picard/MarkDuplicates_RemoveDuplicates-{RemoveDuplicates}_AssumeSorted-{AssumeSorted}/{id}.metrics.txt",
    wildcard_constraints:
        RemoveDuplicates="true|false",
        AssumeSorted="true|false"
    conda:
        "../envs/picard.yaml"
    shell:
        """
        picard\
            MarkDuplicates\
            INPUT={input.bam}\
            OUTPUT={output.bam}\
            METRICS_FILE={output.metrics}\
            REMOVE_DUPLICATES={wildcards.RemoveDuplicates}\
            ASSUME_SORTED={wildcards.AssumeSorted}\
            VALIDATION_STRINGENCY=SILENT
        """


