rule gunzip_extra:
    """
    Created:
        2019-01-30 11:32:31
    Test:
        out/gunzip/_-c/rsync/_-aP/hgdownload.cse.ucsc.edu/goldenPath/hg19/database/tfbsConsSites.txt
        out/gunzip/to-stdout/rsync/_-aP/hgdownload.cse.ucsc.edu/goldenPath/hg19/database/tfbsConsSites.txt
    """
    input:
        "out/{filler}.gz"
    output:
        "out/{tool}{extra}/{filler}"
    log:
        "out/{tool}{extra}/{filler}.log"
    benchmark:
        "out/{tool}{extra}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="gunzip/",
        #extra= Maybe I should add here exclusion if extra start with 'merge' so other rules below are not ambiguous
    shell:
        "gunzip {params.extra} {input} > {output} 2> {log}"

rule gzip_extra:
    """
    Created:
        2017-09-14 15:07:05
    Test:
        "out/gzip/to-stdout/awk/extract_feature_in_repeatMasker/mm9/AT_rich.bed.gz"
    """
    input:
        "out/{filler}"
    output:
        "out/{tool}{extra}/{filler}.gz"
    log:
        "out/{tool}{extra}/{filler}.log"
    benchmark:
        "out/{tool}{extra}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="gzip/",
        #extra= Maybe I should add here exclusion if extra start with 'merge' so other rules below are not ambiguous
    shell:
        "gzip {params.extra} {input} > {output} 2> {log}"

rule gunzip_merge_lanes_nextseq500_paired_end_legacy:
    """
    Created:
        2017-04-14 14:20:20
    Modified:
        2017-12-08 11:10:37 - Changed to legacy because using 'pe' instead of 'paired_end' is more consistent with what is used for alignement and trimming, i.e. shared wildcards between rules.
    Aim:
        Merge the 4 lanes produced by Illumina NextSeq500. 
    Test:
        out/gunzip/merge_lanes_nextseq500_paired_end/run140/H4K5ac-Nut-WT_R1.fastq
    """
    input:
        fastq_lane1_mate1 = "out/{filler}/L1_R1.fastq.gz",
        fastq_lane1_mate2 = "out/{filler}/L1_R2.fastq.gz",
        fastq_lane2_mate1 = "out/{filler}/L2_R1.fastq.gz",
        fastq_lane2_mate2 = "out/{filler}/L2_R2.fastq.gz",
        fastq_lane3_mate1 = "out/{filler}/L3_R1.fastq.gz",
        fastq_lane3_mate2 = "out/{filler}/L3_R2.fastq.gz",
        fastq_lane4_mate1 = "out/{filler}/L4_R1.fastq.gz",
        fastq_lane4_mate2 = "out/{filler}/L4_R2.fastq.gz"
    output:
        fastq_mate1="out/gunzip/merge_lanes_nextseq500_paired_end/{filler}_R1.fastq",
        fastq_mate2="out/gunzip/merge_lanes_nextseq500_paired_end/{filler}_R2.fastq"
    threads:
        1
    shell:
        """
        gunzip --to-stdout {input.fastq_lane1_mate1} > \
            {output.fastq_mate1}

        gunzip --to-stdout {input.fastq_lane1_mate2} > \
            {output.fastq_mate2}

        gunzip --to-stdout {input.fastq_lane2_mate1} >> \
            {output.fastq_mate1}

        gunzip --to-stdout {input.fastq_lane2_mate2} >> \
            {output.fastq_mate2}

        gunzip --to-stdout {input.fastq_lane3_mate1} >> \
            {output.fastq_mate1}

        gunzip --to-stdout {input.fastq_lane3_mate2} >> \
            {output.fastq_mate2}

        gunzip --to-stdout {input.fastq_lane4_mate1} >> \
            {output.fastq_mate1}

        gunzip --to-stdout {input.fastq_lane4_mate2} >> \
            {output.fastq_mate2}
        """

rule gunzip_merge_lanes_nextseq500_pe:
    """
    Created:
        2017-12-08 11:10:01
    Aim:
        Merge the 4 lanes produced by Illumina NextSeq500. 
    Test:
        out/gunzip/merge_lanes_nextseq500_paired_end/run140/H4K5ac-Nut-WT_R1.fastq
    """
    input:
        fastq_lane1_mate1 = "out/{filler}/L1_R1.fastq.gz",
        fastq_lane1_mate2 = "out/{filler}/L1_R2.fastq.gz",
        fastq_lane2_mate1 = "out/{filler}/L2_R1.fastq.gz",
        fastq_lane2_mate2 = "out/{filler}/L2_R2.fastq.gz",
        fastq_lane3_mate1 = "out/{filler}/L3_R1.fastq.gz",
        fastq_lane3_mate2 = "out/{filler}/L3_R2.fastq.gz",
        fastq_lane4_mate1 = "out/{filler}/L4_R1.fastq.gz",
        fastq_lane4_mate2 = "out/{filler}/L4_R2.fastq.gz"
    output:
        fastq_mate1="out/gunzip/merge_lanes_nextseq500_pe/{filler}_R1.fastq",
        fastq_mate2="out/gunzip/merge_lanes_nextseq500_pe/{filler}_R2.fastq"
    threads:
        1
    shell:
        """
        gunzip --to-stdout {input.fastq_lane1_mate1} > \
            {output.fastq_mate1}

        gunzip --to-stdout {input.fastq_lane1_mate2} > \
            {output.fastq_mate2}

        gunzip --to-stdout {input.fastq_lane2_mate1} >> \
            {output.fastq_mate1}

        gunzip --to-stdout {input.fastq_lane2_mate2} >> \
            {output.fastq_mate2}

        gunzip --to-stdout {input.fastq_lane3_mate1} >> \
            {output.fastq_mate1}

        gunzip --to-stdout {input.fastq_lane3_mate2} >> \
            {output.fastq_mate2}

        gunzip --to-stdout {input.fastq_lane4_mate1} >> \
            {output.fastq_mate1}

        gunzip --to-stdout {input.fastq_lane4_mate2} >> \
            {output.fastq_mate2}
        """

rule gunzip_merge_lanes_nextseq500_pe_raw:
    """
    Created:
        2018-01-04 00:24:54
    Aim:
        Merge the 4 lanes produced by Illumina NextSeq500. 
    Test:
        out/gunzip/merge_lanes_nextseq500_paired_end/run140/H4K5ac-Nut-WT_R1.fastq
    """
    input:
        fastq_lane1_mate1 = "out/{filler}_L001_R1_001.fastq.gz",
        fastq_lane1_mate2 = "out/{filler}_L001_R2_001.fastq.gz",
        fastq_lane2_mate1 = "out/{filler}_L002_R1_001.fastq.gz",
        fastq_lane2_mate2 = "out/{filler}_L002_R2_001.fastq.gz",
        fastq_lane3_mate1 = "out/{filler}_L003_R1_001.fastq.gz",
        fastq_lane3_mate2 = "out/{filler}_L003_R2_001.fastq.gz",
        fastq_lane4_mate1 = "out/{filler}_L004_R1_001.fastq.gz",
        fastq_lane4_mate2 = "out/{filler}_L004_R2_001.fastq.gz"
    output:
        fastq_mate1="out/gunzip/merge_lanes_nextseq500_pe_raw/{filler}_1.fastq",
        fastq_mate2="out/gunzip/merge_lanes_nextseq500_pe_raw/{filler}_2.fastq"
    threads:
        1
    shell:
        """
        gunzip --to-stdout {input.fastq_lane1_mate1} > \
            {output.fastq_mate1}

        gunzip --to-stdout {input.fastq_lane1_mate2} > \
            {output.fastq_mate2}

        gunzip --to-stdout {input.fastq_lane2_mate1} >> \
            {output.fastq_mate1}

        gunzip --to-stdout {input.fastq_lane2_mate2} >> \
            {output.fastq_mate2}

        gunzip --to-stdout {input.fastq_lane3_mate1} >> \
            {output.fastq_mate1}

        gunzip --to-stdout {input.fastq_lane3_mate2} >> \
            {output.fastq_mate2}

        gunzip --to-stdout {input.fastq_lane4_mate1} >> \
            {output.fastq_mate1}

        gunzip --to-stdout {input.fastq_lane4_mate2} >> \
            {output.fastq_mate2}
        """

rule gunzip_already_merged_lane_nextseq500_pe:
    input:
        fastq_mate1 = "out/{filler}_R1_001.fastq.gz",
        fastq_mate2 = "out/{filler}_R2_001.fastq.gz"
    output:
        fastq_mate1 = "out/gunzip/already_merged_lane_nextseq500_pe/{filler}_1.fastq",
        fastq_mate2 = "out/gunzip/already_merged_lane_nextseq500_pe/{filler}_2.fastq"
    shell:
        """
        gunzip --to-stdout {input.fastq_mate1} > {output.fastq_mate1}
        gunzip --to-stdout {input.fastq_mate2} > {output.fastq_mate2}
        """

rule gunzip_already_merged_lane_nextseq500_se:
    input:
    output:
    shell:
        "echo 'never occured before'"

rule gunzip_merge_lanes_nextseq500_single_end_legacy:
    """
    Created:
        2017-05-04 14:18:53
    Modified:
        2017-12-07 18:45:44 - Changed to legacy because 'se' is a pattern shared with following rules and can be used for salva sst project more easily.
    Aim:
        Merge the 4 lanes produced by Illumina NextSeq500. 
    Test:
        out/gunzip/merge_lanes_nextseq500_single_end/ln/rename_run107_tgml/Input.fastq
    """
    input:
        fastq_lane1_mate1 = "out/{filler}/L1_R1.fastq.gz",
        fastq_lane2_mate1 = "out/{filler}/L2_R1.fastq.gz",
        fastq_lane3_mate1 = "out/{filler}/L3_R1.fastq.gz",
        fastq_lane4_mate1 = "out/{filler}/L4_R1.fastq.gz",
    output:
        fastq_mate1="out/gunzip/merge_lanes_nextseq500_single_end/{filler}.fastq",
    threads:
        1
    shell:
        """
        gunzip --to-stdout {input.fastq_lane1_mate1} > \
            {output.fastq_mate1}

        gunzip --to-stdout {input.fastq_lane2_mate1} >> \
            {output.fastq_mate1}

        gunzip --to-stdout {input.fastq_lane3_mate1} >> \
            {output.fastq_mate1}

        gunzip --to-stdout {input.fastq_lane4_mate1} >> \
            {output.fastq_mate1}
        """

rule gunzip_merge_lanes_nextseq500_se:
    """
    Created:
        2017-12-07 18:44:32
    Aim:
        Merge the 4 lanes produced by Illumina NextSeq500. 
    Test:
        out/gunzip/merge_lanes_nextseq500_se/ln/rename_run107_tgml/Input.fastq
    """
    input:
        fastq_lane1_mate1 = "out/{filler}/L1_R1.fastq.gz",
        fastq_lane2_mate1 = "out/{filler}/L2_R1.fastq.gz",
        fastq_lane3_mate1 = "out/{filler}/L3_R1.fastq.gz",
        fastq_lane4_mate1 = "out/{filler}/L4_R1.fastq.gz",
    output:
        fastq_mate1="out/gunzip/merge_lanes_nextseq500_se/{filler}.fastq",
    threads:
        1
    shell:
        """
        gunzip --to-stdout {input.fastq_lane1_mate1} > \
            {output.fastq_mate1}

        gunzip --to-stdout {input.fastq_lane2_mate1} >> \
            {output.fastq_mate1}

        gunzip --to-stdout {input.fastq_lane3_mate1} >> \
            {output.fastq_mate1}

        gunzip --to-stdout {input.fastq_lane4_mate1} >> \
            {output.fastq_mate1}
        """

rule gunzip_merge_lanes_nextseq500_se_raw:
    """
    Created:
        2018-01-02 01:11:37
    Aim:
        Merge the 4 lanes produced by Illumina NextSeq500. 
        Testing merging from raw pattern produced by basespace. Doing that may help in the resolution of bad file naming in sst projet.
    Test:
        out/gunzip/merge_lanes_nextseq500_se_raw/inp/fastq/run146/RD_ATAC-seq_Salvatore-Spicuglia-19020/S001387_TH148_149_CD34pos_1aneg_7neg-19019/TH148-149-CD34pos-1aneg-7neg_S1.fastq
    """
    input:
        fastq_lane1_mate1 = "out/{filler}_L001_R1_001.fastq.gz",
        fastq_lane2_mate1 = "out/{filler}_L002_R1_001.fastq.gz",
        fastq_lane3_mate1 = "out/{filler}_L003_R1_001.fastq.gz",
        fastq_lane4_mate1 = "out/{filler}_L004_R1_001.fastq.gz",
    output:
        fastq_mate1="out/gunzip/merge_lanes_nextseq500_se_raw/{filler}.fastq",
    threads:
        1
    shell:
        """
        gunzip --to-stdout {input.fastq_lane1_mate1} > \
            {output.fastq_mate1}

        gunzip --to-stdout {input.fastq_lane2_mate1} >> \
            {output.fastq_mate1}

        gunzip --to-stdout {input.fastq_lane3_mate1} >> \
            {output.fastq_mate1}

        gunzip --to-stdout {input.fastq_lane4_mate1} >> \
            {output.fastq_mate1}
        """

rule gunzip_merge_illumina_fastq_sets:
    """
    Created:
        2018-01-02 01:11:37
    Aim:
        Merge all samples with same name but different sets by Illumina sequencers.
        These sets are variables depending on sample in the Blueprint project.
    Note:
        Illumina FASTQ files use the following naming scheme:
        <sample name>_<barcode sequence>_L<lane>_R<read number>_<set number>.fastq.gz
    Test:
        out/gunzip/merge_illumina_fastq_sets/fastq/blueprint/8580_LC_TH101_H3K27ac/lane2_8580_ACAGTG_L002_R1.fastq
    """
    input:
        fastq = "inp/{filler}_001.fastq.gz",
    output:
        fastq = "out/gunzip/merge_illumina_fastq_sets/{filler}.fastq",
    threads:
        1
    shell:
        """
        DIRNAME=`dirname {input.fastq}`
        
        rm -f {output.fastq}

        #for FASTQ in `ls -d $DIRNAME/**`
        for FASTQ in `ls -d $DIRNAME/*.fastq.gz`
        do
            echo "Extracting $FASTQ"
            gunzip --to-stdout $FASTQ >> {output.fastq}
        done
        """


#rule gunzip_nut_rnaseq_nut_wrong_naming:
#    """
#    Created:
#        2016-10-27 14h36
#    Warning:
#        I wrote this rule having the file "inp/nut_rnaseq/samples_nut_160317.txt in mind. But Florent warned me the 6th november 2016 about the fact there are no replicates and instead Pachytene and Round spermatid samples.
#    """
#    input:
#        fq=expand("inp/nut_rnaseq/Clean/{id}/pe_file_{mate}.fq.gz", id=["30","31","32", "33"], mate=["1","2"])
#    output:
#        fq=expand("out/fastq/nut_rnaseq/{new_id}_{mate}.fastq", new_id=["Nut-KO_rep1","Nut-KO_rep2","Nut-WT_rep1", "Nut-WT_rep2"], mate=["1","2"])
#    shell:"""
#    gunzip --to-stdout inp/nut_rnaseq/Clean/30/pe_file_1.fq.gz > out/fastq/nut_rnaseq/Nut-KO_rep1_1.fastq
#    gunzip --to-stdout inp/nut_rnaseq/Clean/31/pe_file_1.fq.gz > out/fastq/nut_rnaseq/Nut-KO_rep2_1.fastq
#    gunzip --to-stdout inp/nut_rnaseq/Clean/32/pe_file_1.fq.gz > out/fastq/nut_rnaseq/Nut-WT_rep1_1.fastq
#    gunzip --to-stdout inp/nut_rnaseq/Clean/33/pe_file_1.fq.gz > out/fastq/nut_rnaseq/Nut-WT_rep2_1.fastq
#    gunzip --to-stdout inp/nut_rnaseq/Clean/30/pe_file_2.fq.gz > out/fastq/nut_rnaseq/Nut-KO_rep1_2.fastq
#    gunzip --to-stdout inp/nut_rnaseq/Clean/31/pe_file_2.fq.gz > out/fastq/nut_rnaseq/Nut-KO_rep2_2.fastq
#    gunzip --to-stdout inp/nut_rnaseq/Clean/32/pe_file_2.fq.gz > out/fastq/nut_rnaseq/Nut-WT_rep1_2.fastq
#    gunzip --to-stdout inp/nut_rnaseq/Clean/33/pe_file_2.fq.gz > out/fastq/nut_rnaseq/Nut-WT_rep2_2.fastq
#    """

rule gunzip_nut_rnaseq:
    """
    Created:
        2016-10-27 14h36
    Modified:
        2016-11-8 10h30 - Correct naming of files.
        2017-09-25 18:23:44 - output directory is now in gunzip and not in fastq.

    Correct naming:
        label   files               P_R WT_KO 
        30      30_mm10_counts.txt  P   KO
        31      31_mm10_counts.txt  R   KO
        32      32_mm10_counts.txt  P   WT
        33      33_mm10_counts.txt  R   WT

    Warning: I wrote this rule having the file "inp/nut_rnaseq/samples_nut_160317.txt in mind. But Florent warned me the 6th november 2016 about the fact there are no replicates and instead Pachytene and Round spermatid samples.
    """
    input:
        fq=expand("inp/nut_rnaseq/Clean/{id}/pe_file_{mate}.fq.gz", id=["30","31","32", "33"], mate=["1","2"])
    output:
        fq=expand("out/gunzip/nut_rnaseq/{new_id}_{mate}.fastq", new_id=["Nut-P-KO","Nut-R-KO","Nut-P-WT", "Nut-R-WT"], mate=["1","2"])
    shell:
        """
        gunzip --to-stdout inp/nut_rnaseq/Clean/30/pe_file_1.fq.gz > out/gunzip/nut_rnaseq/Nut-P-KO_1.fastq
        gunzip --to-stdout inp/nut_rnaseq/Clean/31/pe_file_1.fq.gz > out/gunzip/nut_rnaseq/Nut-R-KO_1.fastq
        gunzip --to-stdout inp/nut_rnaseq/Clean/32/pe_file_1.fq.gz > out/gunzip/nut_rnaseq/Nut-P-WT_1.fastq
        gunzip --to-stdout inp/nut_rnaseq/Clean/33/pe_file_1.fq.gz > out/gunzip/nut_rnaseq/Nut-R-WT_1.fastq
        gunzip --to-stdout inp/nut_rnaseq/Clean/30/pe_file_2.fq.gz > out/gunzip/nut_rnaseq/Nut-P-KO_2.fastq
        gunzip --to-stdout inp/nut_rnaseq/Clean/31/pe_file_2.fq.gz > out/gunzip/nut_rnaseq/Nut-R-KO_2.fastq
        gunzip --to-stdout inp/nut_rnaseq/Clean/32/pe_file_2.fq.gz > out/gunzip/nut_rnaseq/Nut-P-WT_2.fastq
        gunzip --to-stdout inp/nut_rnaseq/Clean/33/pe_file_2.fq.gz > out/gunzip/nut_rnaseq/Nut-R-WT_2.fastq
        """

rule merge_lanes_nextseq500_paired_end_legacy:
    """
    Created:
        2016
    Modifed:
        2017-03-08 11:30:55 - Updated paths.
        2017-04-14 14:20:40 - Changed to legacy.
    Aim:
        Merge the 4 lanes produced by Illumina NextSeq500
    Test:
        out/merge_lanes/run140/H4K5ac-Nut-WT_R1.fastq
    """
    input:
        fastq_lane1_mate1 = "inp/fastq/{sample}/L1_R1.fastq.gz",
        fastq_lane1_mate2 = "inp/fastq/{sample}/L1_R2.fastq.gz",
        fastq_lane2_mate1 = "inp/fastq/{sample}/L2_R1.fastq.gz",
        fastq_lane2_mate2 = "inp/fastq/{sample}/L2_R2.fastq.gz",
        fastq_lane3_mate1 = "inp/fastq/{sample}/L3_R1.fastq.gz",
        fastq_lane3_mate2 = "inp/fastq/{sample}/L3_R2.fastq.gz",
        fastq_lane4_mate1 = "inp/fastq/{sample}/L4_R1.fastq.gz",
        fastq_lane4_mate2 = "inp/fastq/{sample}/L4_R2.fastq.gz"
    output:
        fastq_mate1="out/merge_lanes/{sample}_R1.fastq",
        fastq_mate2="out/merge_lanes/{sample}_R2.fastq"
    threads:
        1
    shell:
        """
        gunzip --to-stdout {input.fastq_lane1_mate1} > \
            out/merge_lanes/{wildcards.sample}_R1.fastq

        gunzip --to-stdout {input.fastq_lane1_mate2} > \
            out/merge_lanes/{wildcards.sample}_R2.fastq

        gunzip --to-stdout {input.fastq_lane2_mate1} >> \
            out/merge_lanes/{wildcards.sample}_R1.fastq

        gunzip --to-stdout {input.fastq_lane2_mate2} >> \
            out/merge_lanes/{wildcards.sample}_R2.fastq

        gunzip --to-stdout {input.fastq_lane3_mate1} >> \
            out/merge_lanes/{wildcards.sample}_R1.fastq

        gunzip --to-stdout {input.fastq_lane3_mate2} >> \
            out/merge_lanes/{wildcards.sample}_R2.fastq

        gunzip --to-stdout {input.fastq_lane4_mate1} >> \
            out/merge_lanes/{wildcards.sample}_R1.fastq

        gunzip --to-stdout {input.fastq_lane4_mate2} >> \
            out/merge_lanes/{wildcards.sample}_R2.fastq
        """

#rule gunzip:
#    """
#     Explanation for the commented first rule below. It produces sometimes this kind of error:$
#PeriodicWildcardError in line 1 of /gpfs/tagc/home/gcharbonnier/grenoble_project/part2/code/snakemake/rules/gunzip.rules:$
#The value .gz in wildcard id is periodically repeated (out/sickle/pe/-t_sanger_-q_20/sickle/pe_-t_sanger_-q_20/fastq/nut_rnaseq/Nut-WT_rep2_1.fastq.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz.gz). This would lead to an infinite recursion. To avoid this, e.g. restrict the wildcards in this rule to certain values.$
#     """
#    input: "{id}.gz"
#    output: "{id}"
#    shell:"""
#    gunzip --to-stdout {input} > {output}
#    """
#
