rule cat:
    """
    Created:
        2017-03-13 14:56:49
    Aim:
        General purpose concatenation rule.
    Test:
        out/sort/coordinates_bed/cat/mm10_nut_wt_ko_h4k5ac_peaks.bed
        out/sort/coordinates_bed/cat/mm9-P5424-K27ac-WT-KO-IK-sharp-peaks.bed
        out/sort/coordinates_bed/cat/mm10-R-H4K5ac-bu-peaks.bed
        out/sort/coordinates_bed/cat/hg38-macs2-peaks-H3K27ac-thymus.bed
    """
    input:
        lambda wildcards: eval(config['ids'][wildcards.cat_id])
    output:
        "out/cat/{cat_id}"
    shell:
        "cat {input} > {output}"


rule cat_merge_2_lanes_pe:
    """
    Created:
        2021-02-05 16:18:13    
    Aim:
        Merge the 4 lanes produced by Illumina NextSeq500.·
        Testing merging from raw pattern produced by basespace. Doing that may help in the resolution of bad file naming in sst projet.
    Test:
        out/cat/merge_lanes_nextseq500_se/inp/fastq/run146/RD_ATAC-seq_Salvatore-Spicuglia-19020/S001387_TH148_149_CD34pos_1aneg_7neg-19019/TH148-149-CD34pos-1aneg-7neg_S1.fastq.gz
    """
    input:
        fastq_lane1 = "out/{filler}_L001_R{n}_001.fastq.gz",
        fastq_lane2 = "out/{filler}_L002_R{n}_001.fastq.gz",
    output:
        fastq = "out/cat/merge_2_lanes_pe/{filler}_R{n}.fastq.gz"
    wildcard_constraints:
        n = "1|2"
    threads:
        1
    shell:
        "cat {input} > {output}"

rule cat_merge_lanes_nextseq500_se:
    """
    Created:
        2018-01-02 01:11:37
    Aim:
        Merge the 4 lanes produced by Illumina NextSeq500.·
        Testing merging from raw pattern produced by basespace. Doing that may help in the resolution of bad file naming in sst projet.
    Test:
        out/cat/merge_lanes_nextseq500_se/inp/fastq/run146/RD_ATAC-seq_Salvatore-Spicuglia-19020/S001387_TH148_149_CD34pos_1aneg_7neg-19019/TH148-149-CD34pos-1aneg-7neg_S1.fastq.gz
    """
    input:
        fastq_lane1_mate1 = "out/{filler}_L001_R1_001.fastq.gz",
        fastq_lane2_mate1 = "out/{filler}_L002_R1_001.fastq.gz",
        fastq_lane3_mate1 = "out/{filler}_L003_R1_001.fastq.gz",
        fastq_lane4_mate1 = "out/{filler}_L004_R1_001.fastq.gz",
    output:
        fastq_mate1="out/cat/merge_lanes_nextseq500_se/{filler}.fastq.gz",
    threads:
        1
    shell:
        "cat {input} > {output}"

rule cat_merge_lanes_nextseq500_pe_LN_RM:
    """
    Created:
        2017-04-14 14:20:20
    Modified:
        2017-12-08 11:10:37 - Changed to legacy because using 'pe' instead of 'paired_end' is more consistent with what is used for alignement and trimming, i.e. shared wildcards between rules.
    Aim:
        Merge the 4 lanes produced by Illumina NextSeq500.·
    Test:
        out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run187/MNase_Spm_WT_band1_R1.fastq.gz
        out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-KO-Rep1_R1.fastq.gz

samtools/view_bSh/java/select_subpopulations_from_bam_lmin-130_lmax-170/samtools/merge_three_samples/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/run163_run167_run184_run187_run205/s1-MNase_Spm_WT_band1_s2-MNase_Spm_WT_band3_s3-MNase_Spm_WT_band6.bed

    """
    input:
        fastq_mate1 = expand("out/{{filler}}/L{lane}_R1.fastq.gz", lane=["1","2","3","4"]),
        fastq_mate2 = expand("out/{{filler}}/L{lane}_R2.fastq.gz", lane=["1","2","3","4"])
    output:
        fastq_mate1="out/cat/merge_lanes_nextseq500_pe_LN_RM/{filler}_R1.fastq.gz",
        fastq_mate2="out/cat/merge_lanes_nextseq500_pe_LN_RM/{filler}_R2.fastq.gz"
    shell:
        "cat {input.fastq_mate1} > {output.fastq_mate1}; "
        "cat {input.fastq_mate2} > {output.fastq_mate2}"

rule cat_merge_lanes_nextseq500_pe_L00N_RM_001:
    """
    Created:
        2019-09-03 16:02:42
    Aim:
        Merge the 4 lanes produced by Illumina NextSeq500.·
    Test:
        out/cat/merge_lanes_nextseq500_pe_L00N_RM_001/ln/updir/mw-sk/inp/fastq/tgml/run192/S001979_1-61061/1_S1__R1.fastq.gz
    """
    input:
        fastq_mate1 = expand("out/{{filler}}L00{lane}_R1_001.fastq.gz", lane=["1","2","3","4"]),
        fastq_mate2 = expand("out/{{filler}}L00{lane}_R2_001.fastq.gz", lane=["1","2","3","4"])
    output:
        fastq_mate1="out/cat/merge_lanes_nextseq500_pe_L00N_RM_001/{filler}_R1.fastq.gz",
        fastq_mate2="out/cat/merge_lanes_nextseq500_pe_L00N_RM_001/{filler}_R2.fastq.gz"
    shell:
        "cat {input.fastq_mate1} > {output.fastq_mate1}; "
        "cat {input.fastq_mate2} > {output.fastq_mate2}"


rule cat_merge_illumina_fastq_sets:
    """
    Created:
		2020-04-26 12:56:40
    Aim:
        Merge all samples with same name but different sets by Illumina sequencers.
        These sets are variables depending on samples so the rule is looking for the first one and concatenate all those with the same sample prefix.
    Note:
        Illumina FASTQ files use the following naming scheme:
        <sample name>_<barcode sequence>_L<lane>_R<read number>_<set number>.fastq.gz
    Test:
        out/cat/merge_illumina_fastq_sets/ln/updir/mw/inp/fastq/blueprint/8580_LC_TH101_H3K27ac/lane2_8580_ACAGTG_L002_R1.fastq.gz
    """
    input:
        fastq = "out/{filler}_001.fastq.gz",
    output:
        fastq = "out/cat/merge_illumina_fastq_sets/{filler}.fastq.gz",
    threads:
        1
    shell:
        """
        SAMPLE_PREFIX=`echo {input} | sed 's/_001.fastq.gz//'`
        cat $SAMPLE_PREFIX*.fastq.gz > {output}
        """

rule cat_nico:
    """
    Created:
        2017-11-08 11:13:25
    Aim:
        Concatenation rule for Nicolas.
        SampleA_rep1_R1 avec SampleA_rep2_R1
        puis SampleA_rep1_R2 avec SampleA_rep2_R2.
        Et ainsi de suite avec SampleB, C et 250 suivant.
    Use:
        Define SAMPLE variable with a global wildcard:
            SAMPLE, = glob_wildcards('inp/Sample{sample}_rep1_R1.fq')
        Then use this in a snakemake input:
            expand('out/cat_nico/Sample{sample}_R{mate}.fq', sample=SAMPLE, mate=['1','2'])
    """
    input:
        fq_rep1="inp/Sample{sample}_rep1_R{mate}.fq",
        fq_rep2="inp/Sample{sample}_rep2_R{mate}.fq"
    output:
        "out/cat_nico/Sample{sample}_R{mate}.fq"
    shell:
        """
        cat {input} > {output}
        """


