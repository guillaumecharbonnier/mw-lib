def input_igv_batch_dependencies(wildcards):
    """
    Created:
        2017-11-03 11:09:16
    Aim:
        A function parsing the content of an IGV batch script returning the paths of dependencies. It allows to write batch script directly without writing Snakefile first.
    """
    id=wildcards['filler']
    tex="src/igv/"+id+".txt"

    # Meaning of the regex:
    # ^[^%].* means line that does not start with %, i.e. not commented in latex.
    # \\includegraphics.+{ means includegraphic call with options until start of file path.
    # (.*) is the file path extracted.
    # } is the end of the includegraphic call.
    # This regex implies avoiding additionnal {} in the file path.
    #pattern = re.compile(r"^[^%].*\\includegraphics.+{(.*)}")
    # Not tested and not sure about the [ $] part. I want to get the regex until a space (for argument like name=toto) or end of line if no additionnal argument to load.
    pattern = re.compile(r"^load\s([^\s]*)[\s$]")
    paths = []

    with open (tex, "rt") as infile:
        for line in infile:
            m = pattern.search(line)
            if m:
                path = m.group(1)
                paths.append(path)

    return(paths)

rule igv:
    """
    Created:
        2017-05-16 18:13:03
    Aim:
        Execute an IGV-batch file
    Note:
        To run igv in batch mode without X11 display:
        Solution 1:
        https://groups.google.com/forum/#!topic/igv-help/7TFYL4LbJhU
        (Xvfb :10 &) && DISPLAY=:10 java -Xmx750m -jar igv.jar -b batch.file && killall Xvfb
        Solution 2:
        http://thedusseldorfer.blogspot.fr/2013/09/remotely-plotting-with-igv-even-without.html

    Test:
        out/igv/batch_test_1.done
        out/igv/testing_rseg_on_h4k5ac_nut_data.done
        out/igv/batch/testing_rseg_to_call_genes_keeping_nucleosomes_in_spermatozoa.done
        out/igv/batch/bs_call_blueprint.done
    """
    input:
        txt="src/igv/{filler}.txt",
        dep=input_igv_batch_dependencies
    output:
        done=touch("out/igv/{filler}.done")
    conda:
        "../envs/igv.yaml"
    shell:
        "igv -b {input.txt}"

rule get_inputs_for_igv_batch:
    """
    Created:
        2017-06-06 10:22:57
    """
    input:
        "out/deepTools/bamCoverage_binSize-10_minMappingQuality-0_normalizeUsingRPKM_extendReads-150/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run113_run119_run124/MNS-P-WT.bw",
        "out/deepTools/bamCoverage_binSize-10_minMappingQuality-0_normalizeUsingRPKM_extendReads-150/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run113_run119_run124/MNS-R-WT.bw",
        "out/deepTools/bamCoverage_binSize-10_minMappingQuality-0_normalizeUsingRPKM_extendReads-150/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run113_run125_run126/MNS-SC-WT.bw",
        "out/deepTools/bamCoverage_binSize-10_minMappingQuality-0_normalizeUsingRPKM_extendReads-150/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run113_run125_run126/PSK-SC-WT.bw",
        "out/deepTools/bamCoverage_binSize-10_minMappingQuality-0_normalizeUsingRPKM_extendReads-150/samtools/sort/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-130_lmax-170/samtools/sort/samtools/view_bSh/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/run163/MNase_Spm_WT_band1.bw",
        "out/deepTools/bamCoverage_binSize-10_minMappingQuality-0_normalizeUsingRPKM_extendReads-150/samtools/sort/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-30_lmax-100/samtools/sort/samtools/view_bSh/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/run163/MNase_Spm_WT_band2.bw",
        "out/deepTools/bamCoverage_binSize-10_minMappingQuality-0_normalizeUsingRPKM_extendReads-150/samtools/merge/samtools/sort/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-130_lmax-170/samtools/sort/samtools/view_bSh/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/run163_run167/MNase_Spm_WT_band3.bw",
        "out/deepTools/bamCoverage_binSize-10_minMappingQuality-0_normalizeUsingRPKM_extendReads-150/samtools/merge/samtools/sort/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-100_lmax-130/samtools/sort/samtools/view_bSh/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/run163_run167/MNase_Spm_WT_band4.bw",
        "out/deepTools/bamCoverage_binSize-10_minMappingQuality-0_normalizeUsingRPKM_extendReads-150/samtools/merge/samtools/sort/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-30_lmax-100/samtools/sort/samtools/view_bSh/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/run163_run167/MNase_Spm_WT_band5.bw",
        "out/deepTools/bamCoverage_binSize-10_minMappingQuality-0_normalizeUsingRPKM_extendReads-1/samtools/sort/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-130_lmax-170/samtools/sort/samtools/view_bSh/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/run163/MNase_Spm_WT_band1.bw",
        "out/deepTools/bamCoverage_binSize-10_minMappingQuality-0_normalizeUsingRPKM_extendReads-1/samtools/sort/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-30_lmax-100/samtools/sort/samtools/view_bSh/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/run163/MNase_Spm_WT_band2.bw",
        "out/deepTools/bamCoverage_binSize-10_minMappingQuality-0_normalizeUsingRPKM_extendReads-1/samtools/merge/samtools/sort/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-130_lmax-170/samtools/sort/samtools/view_bSh/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/run163_run167/MNase_Spm_WT_band3.bw",
        "out/deepTools/bamCoverage_binSize-10_minMappingQuality-0_normalizeUsingRPKM_extendReads-1/samtools/merge/samtools/sort/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-100_lmax-130/samtools/sort/samtools/view_bSh/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/run163_run167/MNase_Spm_WT_band4.bw",
        "out/deepTools/bamCoverage_binSize-10_minMappingQuality-0_normalizeUsingRPKM_extendReads-1/samtools/merge/samtools/sort/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-30_lmax-100/samtools/sort/samtools/view_bSh/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/run163_run167/MNase_Spm_WT_band5.bw"

