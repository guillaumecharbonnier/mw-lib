rule deepTools_bamCompare_extra:
    """
    Aim:
        Generic rule to produce bw file of log2ratio of two samples.
    Created;
        2019-02-01 02:58:46
    Test:
        out/deepTools/bamCompare_--scaleFactorsMethod_SES_--binSize_10_--extendReads_150/samtools/index/samtools/sort/samtools/view_sam_to_bam/bowtie2/se_mm10/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR1202037_vs_SRR1202038.bw

    """
    input:
        bam1="out/{filler}/{id_bam1}.bam",
        bam2="out/{filler}/{id_bam2}.bam",
        bai1="out/{filler}/{id_bam1}.bam.bai",
        bai2="out/{filler}/{id_bam2}.bam.bai"
    output:
        bw="out/{tool}{extra}/{filler}/{id_bam1}_vs_{id_bam2}.bw"
    params:
        extra = params_extra
    threads:
        MAX_THREADS
    wildcard_constraints:
        tool='deepTools/bamCompare'
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        bamCompare \
            --bamfile1 {input.bam1} \
            --bamfile2 {input.bam2} \
            --numberOfProcessors {threads} \
            {params.extra} \
            --outFileName {output.bw}
        """


# Legacy past this point but I have to adjust workflows to rule above.

rule deepTools_bamCompare_scaleFactorsMethod_ratio_binsize_extendReads:
    """
    Aim:
        Generic rule to produce bw file of log2ratio of two samples.
    Modified:
        2016-11-8 11h04
    Test:
        out/deepTools/bamCompare_scaleFactorsMethod-SES_ratio-log2_binSize-10_extendReads-150/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run140_run141/H4K5ac-Nut-WT_vs_H4K5ac-Nut-KO.bw
        
    Spike-In test:
        out/deepTools/bamCompare_scaleFactorsMethod-readCount_ratio-subtract_binSize-50_extendReads-/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/ln/alias/fastq/run207/H4K5ac-Nut-WT_vs_Input-Nut-WT.bw
        out/deepTools/bamCompare_scaleFactorsMethod-readCount_ratio-subtract_binSize-50_extendReads-/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/ln/alias/fastq/run207/H4K5ac-Nut-KO_vs_Input-Nut-KO.bw
        out/samtools/idxstats/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_dm6/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/ln/alias/fastq/run207/H4K5ac-Nut-WT.tsv
        out/samtools/idxstats/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_dm6/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/ln/alias/fastq/run207/H4K5ac-Nut-KO.tsv
        out/samtools/flagstat/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_dm6/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/ln/alias/fastq/run207/H4K5ac-Nut-WT.tsv
        out/samtools/flagstat/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_dm6/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/ln/alias/fastq/run207/H4K5ac-Nut-KO.tsv

    """
    input:
        bamCompare="opt/miniconda/envs/deeptools/bin/bamCompare",
        bam1="out/{filler}/{id_bam1}.bam",
        bam2="out/{filler}/{id_bam2}.bam",
        bai1="out/{filler}/{id_bam1}.bam.bai",
        bai2="out/{filler}/{id_bam2}.bam.bai"
    output:
        bw="out/deepTools/bamCompare_scaleFactorsMethod-{scaleFactorsMethod}_ratio-{ratio}_binSize-{binSize}_extendReads-{extendReads}/{filler}/{id_bam1}_vs_{id_bam2}.bw"
    threads:
        MAX_THREADS
    wildcard_constraints:
        ratio="log2|ratio|subtract|add|reciprocal_ratio",
        scaleFactorsMethod="readCount|SES",
        binSize="[0-9]+",
        extendReads="|[0-9]+"
    shell:
        """
        {input.bamCompare} \
            --bamfile1 {input.bam1} \
            --bamfile2 {input.bam2} \
            --numberOfProcessors {threads} \
            --binSize {wildcards.binSize} \
            --ratio {wildcards.ratio} \
            --scaleFactorsMethod {wildcards.scaleFactorsMethod} \
            --extendReads {wildcards.extendReads} \
            --outFileName {output.bw}
        """

rule deepTools_bamCompare_scaleFactorsMethod_ratio_binsize_normalizeUsingRPKM_extendReads:
    """
    Aim:
        Generic rule to produce bw file of log2ratio of two samples. I try RPKM here.
    Created:
        2017-10-25 15:10:32
    Test:
    """
    input:
        bamCompare="opt/miniconda/envs/deeptools/bin/bamCompare",
        bam1="out/{filler}/{id_bam1}.bam",
        bam2="out/{filler}/{id_bam2}.bam",
        bai1="out/{filler}/{id_bam1}.bam.bai",
        bai2="out/{filler}/{id_bam2}.bam.bai"
    output:
        bw="out/deepTools/bamCompare_scaleFactorsMethod-{scaleFactorsMethod}_ratio-{ratio}_binSize-{binSize}_normalizeUsingRPKM_extendReads-{extendReads}/{filler}/{id_bam1}_vs_{id_bam2}.bw"
    threads:
        16
    wildcard_constraints:
        ratio="log2|ratio|subtract|add|reciprocal_ratio",
        scaleFactorsMethod="readCount|SES",
        binSize="[0-9]+",
        extendReads="|[0-9]+"
    shell:
        """
        {input.bamCompare} \
            --bamfile1 {input.bam1} \
            --bamfile2 {input.bam2} \
            --numberOfProcessors {threads} \
            --binSize {wildcards.binSize} \
            --ratio {wildcards.ratio} \
            --scaleFactorsMethod {wildcards.scaleFactorsMethod} \
            --normalizeUsingRPKM \
            --extendReads {wildcards.extendReads} \
            --outFileName {output.bw}
        """


rule deepTools_bamCompare_scaleFactorsMethod_ratio_binSize_minMappingQuality:
    """
    Modified:
        2016-11-8 11h04
        2017-04-11 10:41:23 - Changed to legacy because does not fit the new pattern arguments now separated by '-' and ordered by the tool doc order.
    Aim:
        Generic rule to produce bw file of ratio of two samples, or normalized by input signal.
    Usage:
        expand("out/deepTools/bamCompare/{ratio}_{scaleFactorsMethod}_{binSize}/{filler}/{id_bam1}_vs_{id_bam2}.bw",
            ratio="log2",
            scaleFactorsMethod="readCount",
            binSize="10",
            filler="star/pe_mm9/sickle/pe_-t_sanger_-q_20/fastq/nut_rnaseq",
            id_bam1="Nut-P-WT",
            id_bam2="Nut-R-WT")
    Test:
        out/deepTools/bamCompare/ratiolog2_scaleFactorsMethodSES_binSize20_minMappingQuality0/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/merge_lanes/run140_run141/H4K5ac-Nut-WT_vs_Input-Nut-WT.bw
    """
    input:
        bamCompare="opt/miniconda/envs/deeptools/bin/bamCompare",
        bam1="out/{filler}/{id_bam1}.bam",
        bam2="out/{filler}/{id_bam2}.bam",
        bai1="out/{filler}/{id_bam1}.bam.bai",
        bai2="out/{filler}/{id_bam2}.bam.bai"
    output:
        bw="out/deepTools/bamCompare_scaleFactorsMethod-{scaleFactorsMethod}_ratio-{ratio}_binSize-{binSize}_minMappingQuality-{minMappingQuality}/{filler}/{id_bam1}_vs_{id_bam2}.bw"
    threads: 
        MAX_THREADS
    wildcard_constraints:
        ratio="log2|ratio|subtract|add|reciprocal_ratio",
        scaleFactorsMethod="readCount|SES",
        binSize="[0-9]+",
        minMappingQuality="[0-9]+"
    shell:
        """
        {input.bamCompare} \
            --bamfile1 {input.bam1} \
            --bamfile2 {input.bam2} \
            --numberOfProcessors {threads} \
            --binSize {wildcards.binSize} \
            --minMappingQuality {wildcards.minMappingQuality} \
            --ratio {wildcards.ratio} \
            --scaleFactorsMethod {wildcards.scaleFactorsMethod} \
            --outFileName {output.bw}
        """

rule deepTools_bamCompare_scaleFactorsMethod_ratio_binSize_smoothLength_minMappingQuality:
    """
    Created:
        2017-04-12 13:56:08
    Aim:
        Generic rule to produce bw file of ratio of two samples, or normalized by input signal.
    Test:
        out/deepTools/bamCompare/ratiolog2_scaleFactorsMethodSES_binSize20_minMappingQuality0/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/merge_lanes/run140_run141/H4K5ac-Nut-WT_vs_Input-Nut-WT.bw
    """
    input:
        bamCompare="opt/miniconda/envs/deeptools/bin/bamCompare",
        bam1="out/{filler}/{id_bam1}.bam",
        bam2="out/{filler}/{id_bam2}.bam",
        bai1="out/{filler}/{id_bam1}.bam.bai",
        bai2="out/{filler}/{id_bam2}.bam.bai"
    output:
        bw="out/deepTools/bamCompare_scaleFactorsMethod-{scaleFactorsMethod}_ratio-{ratio}_binSize-{binSize}_smoothLength-{smoothLength}_minMappingQuality-{minMappingQuality}/{filler}/{id_bam1}_vs_{id_bam2}.bw"
    threads: 
        MAX_THREADS
    wildcard_constraints:
        ratio="log2|ratio|subtract|add|reciprocal_ratio",
        scaleFactorsMethod="readCount|SES",
        binSize="[0-9]+",
        smoothLength="[0-9]+",
        minMappingQuality="[0-9]+"
    shell:
        """
        {input.bamCompare} \
            --bamfile1 {input.bam1} \
            --bamfile2 {input.bam2} \
            --numberOfProcessors {threads} \
            --binSize {wildcards.binSize} \
            --smoothLength {wildcards.smoothLength} \
            --minMappingQuality {wildcards.minMappingQuality} \
            --ratio {wildcards.ratio} \
            --scaleFactorsMethod {wildcards.scaleFactorsMethod} \
            --outFileName {output.bw}
        """

rule deepTools_bamCompare_scaleFactorsMethod_ratio_binSize_blackListFileName_smoothLength_minMappingQuality:
    """
    Created:
        2017-04-12 13:56:08
    Aim:
        Generic rule to produce bw file of ratio of two samples, or normalized by input signal.
    """
    input:
        bamCompare="opt/miniconda/envs/deeptools/bin/bamCompare",
        blackListFileName="out/bedtools/complement_hg38/sort/coordinates_bed/awk/extract_main_chr/crossmap/bed_hg19_to_hg38/input/bed/salva/Regions_capture_SE_Starr-seq_Alex_HG19_Merged_cellLine_specific_SE.bed",
        bam1="out/{filler}/{id_bam1}.bam",
        bam2="out/{filler}/{id_bam2}.bam",
        bai1="out/{filler}/{id_bam1}.bam.bai",
        bai2="out/{filler}/{id_bam2}.bam.bai"
    output:
        bw="out/deepTools/bamCompare_scaleFactorsMethod-{scaleFactorsMethod}_ratio-{ratio}_binSize-{binSize}_blackListFileName-no-super-enhancers-hg38_smoothLength-{smoothLength}_minMappingQuality-{minMappingQuality}/{filler}/{id_bam1}_vs_{id_bam2}.bw"
    threads:
        MAX_THREADS
    wildcard_constraints:
        ratio="log2|ratio|subtract|add|reciprocal_ratio",
        scaleFactorsMethod="readCount|SES",
        binSize="[0-9]+",
        smoothLength="[0-9]+",
        minMappingQuality="[0-9]+"
    shell:
        """
        {input.bamCompare} \
            --bamfile1 {input.bam1} \
            --bamfile2 {input.bam2} \
            --numberOfProcessors {threads} \
            --binSize {wildcards.binSize} \
            --blackListFileName {input.blackListFileName} \
            --smoothLength {wildcards.smoothLength} \
            --minMappingQuality {wildcards.minMappingQuality} \
            --ratio {wildcards.ratio} \
            --scaleFactorsMethod {wildcards.scaleFactorsMethod} \
            --outFileName {output.bw}
        """

rule deepTools_bamCompare_scaleFactorsMethod_ratio_binSize_smoothLength_extendReads:
    """
    Created:
        2017-11-03 11:05:45
    Aim:
        Generic rule to produce bw file of ratio of two samples, or normalized by input signal.
    Test:
        out/deepTools/bamCompare_scaleFactorsMethod-readCount_ratio-log2_binSize-40_smoothLength-80_extendReads-/samtools/merge/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_mm9/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run140_run141/H4K5ac-Nut-WT_vs_H4K5ac-Nut-KO.bw
    """
    input:
        bamCompare="opt/miniconda/envs/deeptools/bin/bamCompare",
        bam1="out/{filler}/{id_bam1}.bam",
        bam2="out/{filler}/{id_bam2}.bam",
        bai1="out/{filler}/{id_bam1}.bam.bai",
        bai2="out/{filler}/{id_bam2}.bam.bai"
    output:
        bw="out/deepTools/bamCompare_scaleFactorsMethod-{scaleFactorsMethod}_ratio-{ratio}_binSize-{binSize}_smoothLength-{smoothLength}_extendReads-{extendReads}/{filler}/{id_bam1}_vs_{id_bam2}.bw"
    threads: MAX_THREADS
    wildcard_constraints:
        ratio="log2|ratio|subtract|add|reciprocal_ratio",
        scaleFactorsMethod="readCount|SES",
        binSize="[0-9]+",
        smoothLength="[0-9]+",
        extendReads="|[0-9]+"
    shell:
        """
        {input.bamCompare} \
            --bamfile1 {input.bam1} \
            --bamfile2 {input.bam2} \
            --numberOfProcessors {threads} \
            --binSize {wildcards.binSize} \
            --smoothLength {wildcards.smoothLength} \
            --ratio {wildcards.ratio} \
            --scaleFactorsMethod {wildcards.scaleFactorsMethod} \
            --extendReads {wildcards.extendReads} \
            --outFileName {output.bw}
        """

rule deepTools_bamCompare_scaleFactorsMethod_ratio_binSize_smoothLength_minMappingQuality_extendReads:
    """
    Created:
        2017-05-18 15:56:10
    Aim:
        Generic rule to produce bw file of ratio of two samples, or normalized by input signal.
    Test:
    """
    input:
        bamCompare="opt/miniconda/envs/deeptools/bin/bamCompare",
        bam1="out/{filler}/{id_bam1}.bam",
        bam2="out/{filler}/{id_bam2}.bam",
        bai1="out/{filler}/{id_bam1}.bam.bai",
        bai2="out/{filler}/{id_bam2}.bam.bai"
    output:
        bw="out/deepTools/bamCompare_scaleFactorsMethod-{scaleFactorsMethod}_ratio-{ratio}_binSize-{binSize}_smoothLength-{smoothLength}_minMappingQuality-{minMappingQuality}_extendReads-{extendReads}/{filler}/{id_bam1}_vs_{id_bam2}.bw"
    threads: 16
    wildcard_constraints:
        ratio="log2|ratio|subtract|add|reciprocal_ratio",
        scaleFactorsMethod="readCount|SES",
        binSize="[0-9]+",
        smoothLength="[0-9]+",
        minMappingQuality="[0-9]+",
        extendReads="|[0-9]+"
    shell:
        """
        {input.bamCompare} \
            --bamfile1 {input.bam1} \
            --bamfile2 {input.bam2} \
            --numberOfProcessors {threads} \
            --binSize {wildcards.binSize} \
            --smoothLength {wildcards.smoothLength} \
            --minMappingQuality {wildcards.minMappingQuality} \
            --ratio {wildcards.ratio} \
            --scaleFactorsMethod {wildcards.scaleFactorsMethod} \
            --extendReads {wildcards.extendReads} \
            --outFileName {output.bw}
        """

rule deepTools_bamCompare_ratio_scaleFactorsMethod_binsize_legacy:
    """
    Generic rule to produce bw file of log2ratio of two samples.
    
    Modified: 2016-11-8 11h04

    Usage:
        expand("out/deepTools/bamCompare/{ratio}_{scaleFactorsMethod}_{binSize}/{filler}/{id_bam1}_vs_{id_bam2}.bw",
            ratio="log2",
            scaleFactorsMethod="readCount",
            binSize="10",
            filler="star/pe_mm9/sickle/pe_-t_sanger_-q_20/fastq/nut_rnaseq",
            id_bam1="Nut-P-WT",
            id_bam2="Nut-R-WT")
    Test:
        out/deepTools/bamCompare/log2_readCount_20/{filler}/{id_bam1}_vs_{id_bam2}.bw
    """
    input:
        bamCompare="opt/miniconda/envs/py27/bin/bamCompare",
        bam1="out/{filler}/{id_bam1}.bam",
        bam2="out/{filler}/{id_bam2}.bam"
    output:
        bw="out/deepTools/bamCompare/{ratio}_{scaleFactorsMethod}_{binSize}/{filler}/{id_bam1}_vs_{id_bam2}.bw"
    threads: 16
    wildcard_constraints:
        ratio="log2|ratio|subtract|add|reciprocal_ratio",
        scaleFactorsMethod="readCount|SES",
        binSize="[0-9]+"
    shell:
        """
        {input.bamCompare} \
            --bamfile1 {input.bam1} \
            --bamfile2 {input.bam2} \
            --numberOfProcessors {threads} \
            --binSize {wildcards.binSize} \
            --ratio {wildcards.ratio} \
            --scaleFactorsMethod {wildcards.scaleFactorsMethod} \
            --outFileName {output.bw}
        """

rule deepTools_bamCompare_ratio_scaleFactorsMethod_binsize_minMappingQuality_legacy:
    """
    Modified:
        2016-11-8 11h04
        2017-04-11 10:41:23 - Changed to legacy because does not fit the new pattern arguments now separated by '-' and ordered by the tool doc order.
    Aim:
        Generic rule to produce bw file of ratio of two samples, or normalized by input signal.
    Usage:
        expand("out/deepTools/bamCompare/{ratio}_{scaleFactorsMethod}_{binSize}/{filler}/{id_bam1}_vs_{id_bam2}.bw",
            ratio="log2",
            scaleFactorsMethod="readCount",
            binSize="10",
            filler="star/pe_mm9/sickle/pe_-t_sanger_-q_20/fastq/nut_rnaseq",
            id_bam1="Nut-P-WT",
            id_bam2="Nut-R-WT")
    Test:
        out/deepTools/bamCompare/ratiolog2_scaleFactorsMethodSES_binSize20_minMappingQuality0/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/merge_lanes/run140_run141/H4K5ac-Nut-WT_vs_Input-Nut-WT.bw
    """
    input:
        bamCompare="opt/miniconda/envs/deeptools/bin/bamCompare",
        bam1="out/{filler}/{id_bam1}.bam",
        bam2="out/{filler}/{id_bam2}.bam",
        bai1="out/{filler}/{id_bam1}.bam.bai",
        bai2="out/{filler}/{id_bam2}.bam.bai"
    output:
        bw="out/deepTools/bamCompare/ratio{ratio}_scaleFactorsMethod{scaleFactorsMethod}_binSize{binSize}_minMappingQuality{minMappingQuality}/{filler}/{id_bam1}_vs_{id_bam2}.bw"
    threads: 16
    wildcard_constraints:
        ratio="log2|ratio|subtract|add|reciprocal_ratio",
        scaleFactorsMethod="readCount|SES",
        binSize="[0-9]+",
        minMappingQuality="[0-9]+"
    shell:
        """
        {input.bamCompare} \
            --bamfile1 {input.bam1} \
            --bamfile2 {input.bam2} \
            --numberOfProcessors {threads} \
            --binSize {wildcards.binSize} \
            --minMappingQuality {wildcards.minMappingQuality} \
            --ratio {wildcards.ratio} \
            --scaleFactorsMethod {wildcards.scaleFactorsMethod} \
            --outFileName {output.bw}
        """

