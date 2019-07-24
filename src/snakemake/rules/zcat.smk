rule zcat_merge_lanes_nextseq500_se:
    """
    Created:
        2018-01-02 01:11:37
    Aim:
        Merge the 4 lanes produced by Illumina NextSeq500.·
        Testing merging from raw pattern produced by basespace. Doing that may help in the resolution of bad file naming in sst projet.
    Test:
        out/zcat/merge_lanes_nextseq500_se/inp/fastq/run146/RD_ATAC-seq_Salvatore-Spicuglia-19020/S001387_TH148_149_CD34pos_1aneg_7neg-19019/TH148-149-CD34pos-1aneg-7neg_S1.fastq.gz
    """
    input:
        fastq_lane1_mate1 = "out/{filler}_L001_R1_001.fastq.gz",
        fastq_lane2_mate1 = "out/{filler}_L002_R1_001.fastq.gz",
        fastq_lane3_mate1 = "out/{filler}_L003_R1_001.fastq.gz",
        fastq_lane4_mate1 = "out/{filler}_L004_R1_001.fastq.gz",
    output:
        fastq_mate1="out/zcat/merge_lanes_nextseq500_se/{filler}.fastq",
    threads:
        1
    conda:
        "../envs/pigz.yaml"
    shell:
        "zcat {input} > {output}"

