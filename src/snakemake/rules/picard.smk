rule picard_CollectInsertSizeMetrics:
    """
    Created:
        2017-05-06 23:16:37
    Test:
        out/picard/CollectInsertSizeMetrics/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run113_run119_run124/MNS-R-WT.pdf

        p_r_paths=expand("out/picard/CollectInsertSizeMetrics/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run113_run119_run124/{id}.bw", id=["MNS-P-WT","MNS-R-WT"])
        cs_paths=expand( "out/picard/CollectInsertSizeMetrics/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run113_run125_run126/{id}.bw", id=["MNS-SC-WT","PSK-SC-WT"])
        spz_nucleosome_pop_paths=expand("out/picard/CollectInsertSizeMetrics/samtools/sam_to_bam/java/select_subpopulations_from_bam/130_170/samtools/sam_to_bam/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/merge_lanes/run163/MNase_Spm_WT_band{int}.bw", int=["1"])
        spz_nucleosome_pop_paths_merge=expand("out/picard/CollectInsertSizeMetrics/samtools/merge/samtools/sam_to_bam/java/select_subpopulations_from_bam/130_170/samtools/sam_to_bam/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/merge_lanes/run163_run167/MNase_Spm_WT_band{int}.bw", int=["1"])
        spz_small_pop_paths=expand("out/picard/CollectInsertSizeMetrics/samtools/sam_to_bam/java/select_subpopulations_from_bam/30_100/samtools/sam_to_bam/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/merge_lanes/run163/MNase_Spm_WT_band{int}.bw", int=["2"])
        spz_small_pop_paths_merge=expand("out/picard/CollectInsertSizeMetrics/samtools/merge/samtools/sam_to_bam/java/select_subpopulations_from_bam/30_100/samtools/sam_to_bam/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/merge_lanes/run163_run167/MNase_Spm_WT_band{int}.bw", int=["2"])

    """
    input:
        bam="out/{filler}.bam"
    output:
        pdf="out/picard/CollectInsertSizeMetrics/{filler}.pdf",
        txt="out/picard/CollectInsertSizeMetrics/{filler}.txt"
    threads: 1
    conda:
        "../envs/picard.yaml"
    shell:
        """
        picard CollectInsertSizeMetrics I={input.bam} O={output.txt} H={output.pdf}
        """

rule tmp_test_picard_on_mnase_data:
    """
    Created:
        2017-06-06 15:27:26
    """
    input:
        "out/picard/CollectInsertSizeMetrics/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run113_run119_run124/MNS-P-WT.pdf",
        "out/picard/CollectInsertSizeMetrics/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run113_run119_run124/MNS-R-WT.pdf",
        "out/picard/CollectInsertSizeMetrics/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run113_run125_run126/MNS-SC-WT.pdf",
        "out/picard/CollectInsertSizeMetrics/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run113_run125_run126/PSK-SC-WT.pdf",
        "out/picard/CollectInsertSizeMetrics/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band1.pdf",
        "out/picard/CollectInsertSizeMetrics/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band2.pdf",
        "out/picard/CollectInsertSizeMetrics/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band3.pdf",
        "out/picard/CollectInsertSizeMetrics/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band4.pdf",



###### 
# Below are legacy rules.
#####

#rule picard_CollectInsertSizeMetrics_legacy:
#    """
#    Created: 2016-11-22 15h59 - Add new filler pattern
#    """
#    input:
#        bam="out/{filler}.bam",
#        java="opt/jre1.8.0_73/bin/java",
#        picard="opt/picard-tools-2.1.1/picard.jar"
#    output:
#        pdf="result/picard/CollectInsertSizeMetrics/{filler}.pdf",
#        txt="result/picard/CollectInsertSizeMetrics/{filler}.txt"
#    threads: 1
#    shell: """
#    {input.java} -jar {input.picard} CollectInsertSizeMetrics I={input.bam} O={output.txt} H={output.pdf}
#    """
#
#rule picard_CollectInsertSizeMetrics_whole_genome:
#    """
#    Created: 2016-03-21 10h49
#    """
#    input:
#        bam="out/bam/{index}/{exp}/{sample}.bam",
#        java="opt/jre1.8.0_73/bin/java",
#        picard="opt/picard-tools-2.1.1/picard.jar"
#    output:
#        pdf="result/picard/CollectInsertSizeMetrics/whole_genome/{index}/{exp}/insert_size_histogram_{sample}.pdf",
#        txt="result/picard/CollectInsertSizeMetrics/whole_genome/{index}/{exp}/insert_size_metrics_{sample}.txt"
#    threads: 1
#    shell: """
#    {input.java} -jar {input.picard} CollectInsertSizeMetrics I={input.bam} O={output.txt} H={output.pdf}
#    """
#
#rule picard_CollectInsertSizeMetrics_for_feature:
#    """
#    Created: 2016-03-22 17h12
#    """
#    input:
#        bam="out/bedtools/intersect/bam_in_feature/{index}/{exp}/{sample}/{feature}.bam",
#        java="opt/jre1.8.0_73/bin/java",
#        picard="opt/picard-tools-2.1.1/picard.jar"
#    output:
#        pdf="result/picard/CollectInsertSizeMetrics/bam_in_feature/{index}/{exp}/{sample}/insert_size_histogram_{feature}.pdf",
#        txt="result/picard/CollectInsertSizeMetrics/bam_in_feature/{index}/{exp}/{sample}/insert_size_metrics_{feature}.txt"
#    threads: 1
#    shell: """
#    {input.java} -jar {input.picard} CollectInsertSizeMetrics I={input.bam} O={output.txt} H={output.pdf}
#    """
#
#rule picard_CollectInsertSizeMetrics_for_filterQ30:
#    """
#    Created: 2016-03-22 17h12
#    """
#    input:
#        bam="out/samtools/filterQ30/bedtools/intersect/bam_in_feature/{index}/{exp}/{sample}/{feature}_Q30.bam",
#        java="opt/jre1.8.0_73/bin/java",
#        picard="opt/picard-tools-2.1.1/picard.jar"
#    output:
#        pdf="result/picard/CollectInsertSizeMetrics/filterQ30/bam_in_feature/{index}/{exp}/{sample}/insert_size_histogram_{feature}.pdf",
#        txt="result/picard/CollectInsertSizeMetrics/filterQ30/bam_in_feature/{index}/{exp}/{sample}/insert_size_metrics_{feature}.txt"
#    threads: 1
#    shell: """
#    {input.java} -jar {input.picard} CollectInsertSizeMetrics I={input.bam} O={output.txt} H={output.pdf}
#    """
#
#
