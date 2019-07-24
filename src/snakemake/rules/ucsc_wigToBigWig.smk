rule ucsc_wigToBigWig_extra:
    """
    Created:
        2018-01-15 14:47:37
    Aim:
        Converting wig from danpos wiq to bigwig.
    Note:
        -clip is nearly mandatory to solve some issues.
    Test:
        out/ucsc/wigToBigWig_-clip_chrominfo-mm10/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/ln/alias/fastq/run207/pasha/result_PE/AllReads/WIGfs_H3K9ac-Nut-KO_mergedReads_elPairsAndEst143_AThr2_bin50.bw
        out/ucsc/wigToBigWig_-clip_chrominfo-mm10/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/ln/alias/fastq/run207/pasha/result_PE/AllReads/WIGfs_H3K9ac-Nut-WT_mergedReads_elPairsAndEst143_AThr2_bin50.bw
        out/ucsc/wigToBigWig_-clip_chrominfo-mm10/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/ln/alias/fastq/run207/pasha/result_PE/AllReads/WIGfs_H4K5ac-Nut-KO_mergedReads_elPairsAndEst152_AThr3_bin50.bw
        out/ucsc/wigToBigWig_-clip_chrominfo-mm10/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/ln/alias/fastq/run207/pasha/result_PE/AllReads/WIGfs_H4K5ac-Nut-WT_mergedReads_elPairsAndEst153_AThr3_bin50.bw
        out/ucsc/wigToBigWig_-clip_chrominfo-mm10/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/ln/alias/fastq/run207/pasha/result_PE/AllReads/WIGfs_Input-Nut-KO_mergedReads_elPairsAndEst153_AThr3_bin50.bw
        out/ucsc/wigToBigWig_-clip_chrominfo-mm10/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/ln/alias/fastq/run207/pasha/result_PE/AllReads/WIGfs_Input-Nut-WT_mergedReads_elPairsAndEst152_AThr4_bin50.bw
        out/ucsc/wigToBigWig_-clip_chrominfo-hg19-main-chr/gunzip/to-stdout/wget/https/www.genboree.org/EdaccData/Current-Release/experiment-sample/Bisulfite-Seq/Mobilized_CD34_Primary_Cells/BI.Mobilized_CD34_Primary_Cells.Bisulfite-Seq.RO_01549.bw
    """
    input:
        wig = "out/{filler}.wig",
        chrominfo = lambda wildcards: config['ids'][wildcards.chrominfo_id]
    output:
        bw = "out/{tool}{extra}_{chrominfo_id}/{filler}.bw"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="ucsc/wigToBigWig"
    priority:
        3
    conda:
        "../envs/ucsc.yaml"
    shell:
        "wigToBigWig {input.wig} {input.chrominfo} {output.bw} {params.extra}"
