rule samtools_merge:
    """
    Created:
        2018-03-16 18:35:54
    Aim:
        This rule use Samtools merge samples as define in salva_runs.tsv by sample_merge_list value.
    Test:
        out/samtools/merge/bam-hg38-Casero2016-thy3-thy4.bam
        out/samtools/merge/bam-hg19-CapStarr_K562_IFN_rep1-merged.bam
    """
    input:
        bam = lambda wildcards: eval(mwconf['ids'][wildcards.bam_list_id])
    output:
        bam = "out/samtools/merge/{bam_list_id}.bam"
    wildcard_constraints:
        bam_list_id = "bam-[a-zA-Z0-9-]+"
    threads:
        1
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools merge {output.bam} {input.bam}"

rule samtools_merge_two_runs:
    """
    Modified:
        2016-12-1 16h19 - Adapted rule to new patterns.
    Aim:
        This rule use Samtools merge to merge two bams sharing the same filename.
    Test: 
        out/samtools/merge_two_runs/samtools/view_sam_to_bam/bowtie2/pe_mm10/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run187_run205/MNase_Spm_WT_band1.bam
    """
    input:
        bam1="out/{filler}/{run1}/{sample}.bam",
        bam2="out/{filler}/{run2}/{sample}.bam"
    output:
        bam="out/samtools/merge_two_runs/{filler}/{run1}_{run2}/{sample}.bam"
    conda:
        "../envs/samtools.yaml"
    wildcard_constraints:
        run1="run[0-9]+",
        run2="run[0-9]+"
    threads:
        1
    shell:
        "samtools merge {output.bam} {input}"

rule samtools_merge_three_runs:
    """
    Modified:
        2016-12-1 16h19 - Adapted rule to new patterns.
    Aim:
        This rule use Samtools merge to merge three bams sharing the same filename.
    Test:
        out/samtools/merge_three_runs/samtools/view_sam_to_bam/bowtie2/pe_mm10/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run184_run187_run205/MNase_Spm_WT_band1.bam
    """
    input:
        bam1="out/{filler}/{run1}/{sample}.bam",
        bam2="out/{filler}/{run2}/{sample}.bam",
        bam3="out/{filler}/{run3}/{sample}.bam"
    output:
        bam="out/samtools/merge_three_runs/{filler}/{run1}_{run2}_{run3}/{sample}.bam"
    conda:
        "../envs/samtools.yaml"
    wildcard_constraints:
        run1="run[0-9]+",
        run2="run[0-9]+",
        run3="run[0-9]+"
    threads:
        1
    shell:
        "samtools merge {output.bam} {input}"

rule samtools_merge_four_runs:
    """
    Created:
        2017-07-26 15:12:36
    Aim:
        This rule use Samtools merge to merge four bams sharing the same filename.
    Test:
        out/samtools/merge_four_runs/samtools/view_sam_to_bam/bowtie2/pe_mm10/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run163_run184_run187_run205/MNase_Spm_WT_band1.bam
    """
    input:
        bam1="out/{filler}/{run1}/{sample}.bam",
        bam2="out/{filler}/{run2}/{sample}.bam",
        bam3="out/{filler}/{run3}/{sample}.bam",
        bam4="out/{filler}/{run4}/{sample}.bam"
    output:
        bam="out/samtools/merge_four_runs/{filler}/{run1}_{run2}_{run3}_{run4}/{sample}.bam"
    conda:
        "../envs/samtools.yaml"
    wildcard_constraints:
        run1="run[0-9]+",
        run2="run[0-9]+",
        run3="run[0-9]+",
        run4="run[0-9]+"
    threads:
        1
    shell:
        "samtools merge {output.bam} {input}"

rule samtools_merge_five_runs:
    """
    Created:
        2017-07-26 15:12:36
    Aim:
        This rule use Samtools merge to merge five bams sharing the same filename.
    Test:
        out/samtools/merge_five_runs/samtools/view_sam_to_bam/bowtie2/pe_mm10/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run163_run167_run184_run187_run205/MNase_Spm_WT_band1.bam
    """
    input:
        bam1="out/{filler}/{run1}/{sample}.bam",
        bam2="out/{filler}/{run2}/{sample}.bam",
        bam3="out/{filler}/{run3}/{sample}.bam",
        bam4="out/{filler}/{run4}/{sample}.bam",
        bam5="out/{filler}/{run5}/{sample}.bam"
    output:
        bam="out/samtools/merge_five_runs/{filler}/{run1}_{run2}_{run3}_{run4}_{run5}/{sample}.bam"
    conda:
        "../envs/samtools.yaml"
    wildcard_constraints:
        run1="run[0-9]+",
        run2="run[0-9]+",
        run3="run[0-9]+",
        run4="run[0-9]+",
        run5="run[0-9]+"
    threads:
        1
    shell:
        "samtools merge {output.bam} {input}"


rule samtools_merge_three_bam_no_exp:
    """
    Modified:
        2017-02-17 10:34:47 - Adapted rule to new patterns.
    Aim:
        This rule use Samtools merge to merge bams.
        I do not see a smart way to specify a variable number of runs to merge without relying on a tedious input function or multiple rules for 2, 3, 4, etc runs...
    Test:
        "out/samtools/merge/{filler}/{sample1}_{sample2}_{sample3}.bam"
    """
    input:
        bam1="out/{filler}/{sample1}.bam",
        bam2="out/{filler}/{sample2}.bam",
        bam3="out/{filler}/{sample3}.bam"
    output:
        bam="out/samtools/merge/{filler}/{sample1}_{sample2}_{sample3}.bam"
    conda:
        "../envs/samtools.yaml"
    wildcard_constraints:
        sample1="SRR[0-9]+",
        sample2="SRR[0-9]+",
        sample3="SRR[0-9]+"
    threads:
        1
    shell:
        "samtools merge {output.bam} {input}"

rule samtools_merge_two_samples:
    """
    Created:
        2017-08-29 14:53:06
    Aim:
        This rule use Samtools merge to merge bams.
        The difference with the previous rule is that here the bam are merged between samples with different names whereas the previous rules merge bam with same sample names and different runs.
    Test:
        out/samtools/merge_two_samples/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167_run184_run187/s1-MNase_Spm_WT_band1_s2-MNase_Spm_WT_band3.bam

        out/samtools/merge_two_samples/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run140_run141/s1-Input-Nut-WT_s2-Input-Nut-KO.bam
    """
    input:
        bam1="out/{filler}/{sample1}.bam",
        bam2="out/{filler}/{sample2}.bam"
    output:
        bam="out/samtools/merge_two_samples/{filler}/s1-{sample1}_s2-{sample2}.bam"
    conda:
        "../envs/samtools.yaml"
    threads:
        1
    shell:
        "samtools merge {output.bam} {input}"

rule samtools_merge_three_samples:
    """
    Created:
        2017-08-29 15:15:38
    Aim:
        This rule use Samtools merge to merge bams.
        The difference with the previous rule is that here the bam are merged between samples with different names whereas the previous rules merge bam with same sample names and different runs.
    Test:
        out/samtools/merge_three_samples/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167_run184_run187/s1-MNase_Spm_WT_band1_s2-MNase_Spm_WT_band3_s3-MNase_Spm_WT_band6.bam
        out/samtools/merge_three_samples/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167_run184/s1-MNase_Spm_WT_band2_s2-MNase_Spm_WT_band4_s3-MNase_Spm_WT_band5.bam
        out/samtools/merge_three_samples/inp/bam/hg38/H3K27ac/thymus/Blueprint/s1-Th118_LC_Input_s2-Th91_LC_input_s3-Th91_SP4_input.bam
        out/samtools/merge_three_samples/samtools/merge_five_runs/samtools/view_sam_to_bam/bowtie2/pe_mm10/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run163_run167_run184_run187_run205/s1-MNase_Spm_WT_band1_s2-MNase_Spm_WT_band3_s3-MNase_Spm_WT_band6.bam
        out/samtools/sort/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-130_lmax-170/out/samtools/merge_three_samples/samtools/merge_five_runs/samtools/view_sam_to_bam/bowtie2/pe_mm10/ln/pe_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run163_run167_run184_run187_run205/s1-MNase_Spm_WT_band1_s2-MNase_Spm_WT_band3_s3-MNase_Spm_WT_band6.bam out/samtools/sort/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-30_lmax-100/samtools/merge_three_samples/samtools/merge_three_runs/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/pe_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/run163_run167_run184/s1-MNase_Spm_WT_band2_s2-MNase_Spm_WT_band4_s3-MNase_Spm_WT_band5.bed

    """
    input:
        bam1="out/{filler}/{sample1}.bam",
        bam2="out/{filler}/{sample2}.bam",
        bam3="out/{filler}/{sample3}.bam"
    output:
        bam="out/samtools/merge_three_samples/{filler}/s1-{sample1}_s2-{sample2}_s3-{sample3}.bam"
    conda:
        "../envs/samtools.yaml"
    threads:
        1
    shell:
        "samtools merge {output.bam} {input}"

rule samtools_merge_four_samples:
    """
    Created:
        2017-09-13 09:06:21
    Aim:
        This rule use Samtools merge to merge four bams.
        The difference with the previous rule is that here the bam are merged between samples with different names whereas the previous rules merge bam with same sample names and different runs.
    """
    input:
        bam1="out/{filler}/{sample1}.bam",
        bam2="out/{filler}/{sample2}.bam",
        bam3="out/{filler}/{sample3}.bam",
        bam4="out/{filler}/{sample4}.bam"
    output:
        bam="out/samtools/merge_four_samples/{filler}/s1-{sample1}_s2-{sample2}_s3-{sample3}_s4-{sample4}.bam"
    conda:
        "../envs/samtools.yaml"
    threads:
        1
    shell:
        "samtools merge {output.bam} {input}"

