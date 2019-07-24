"""
Created: 2016-10-14 10h17 - This file will contain multiple ways to compute FRIP-like metric to answer Saadi's main question on H4K5ac Nut WT/KO ChIP-seq data.
"""

rule chipqc:
    """
    Created:
        2019-03-03 19:47:03
    Aim:
        Test this QC tool.
    Note:
        Conda fail to install the environment:
        ERROR conda.core.link:_execute(502): An error occurred while installing package 'bioconda::bioconductor-genomeinfodbdata-1.2.0-r351_0'.
        LinkError: post-link script failed for package bioconda::bioconductor-genomeinfodbdata-1.2.0-r351_0
    Test:
        out/chipqc/sample/samtools/sort/samtools/view_bSh/bowtie2/se_mm10/sickle/se_-t_sanger_-q_20/sra-tools/fastqdump_se/wget/ftp_trace_ncbi_sra/SRR3126243.bam
    """
    input:
        bam="out/{filler}.bam"
    output:
        pdf="out/chipqc/sample/{filler}.pdf"
    conda:
        "../envs/chipqc.yaml"
    script:
        "src/r/script/chipqc_sample.R"

rule frip_tmp:
    """
    Note: when subsampling with samtools view, paired-end information is not conserved, thus I require "BAM" instead of "BAMPE"
    """
    input:
        pdf="result/frip/from_macs2_output/macs2/callpeak_broad_BAMPE_mm_SPMRFalse/samtools/view/subsampling_bam/seed1_55000000pf/picard/MarkDuplicates/bam/mm9/merge_run140_run141/H4K5ac-Nut-WT_over_Input-Nut-WT_VS_H4K5ac-Nut-KO_over_Input-Nut-KO.pdf",
        control="result/frip/from_macs2_output/macs2/callpeak_broad_BAMPE_mm_SPMRFalse/samtools/view/subsampling_bam/seed7_30000000pf/bam/mm9/merge_run140_run141/H4K5ac-Nut-WT_over_Input-Nut-WT_VS_Input-Nut-KO_over_Input-Nut-WT.pdf"

rule r_peak_profiles_from_macs2_output:
    """
    Created:
        2017-10-05 16:03:54
    Aim:
        Produce plots of length and pileup to look the profile of peaks produced by MACS2.
    Test:
        out/r/peak_profile_from_macs2_output/macs2/callpeak_broad_format-BAMPE_gsize-mm_SPMR/samtools/view/subsampling_bam/seed7_30000000pf/picard/MarkDuplicates/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run140_run141/H4K5ac-Nut-WT_over_Input-Nut-WT.pdf

    """
    input:
        xls="out/{filler}_peaks.xls",
        code="src/r/script/peak_profile_from_macs2_output.R",
        Rscript="opt/miniconda/envs/ggrepel/bin/Rscript"
    output:
        pdf="out/r/peak_profile_from_macs2_output/{filler}.pdf"
    shell:
        """
        {input.Rscript} {input.code} -i {input.xls} -o {output.pdf}
        """

rule r_frip_from_macs2_output:
    """
    Created: 2016-11-14 9h56

    xls_wt="out/macs2/callpeak_broad_BAMPE_mm_SPMRTrue/bam/mm9/merge_run140_run141/H4K5ac-Nut-WT_over_Input-Nut-WT_peaks.xls",
    xls_ko="out/macs2/callpeak_broad_BAMPE_mm_SPMRTrue/bam/mm9/merge_run140_run141/H4K5ac-Nut-KO_over_Input-Nut-KO_peaks.xls"

    Test:
        out/r/frip_from_macs2_output/macs2/callpeak_broad_format-BAMPE_gsize-mm_SPMR/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run140_run141/H4K5ac-Nut-WT_over_Input-Nut-WT_VS_Input-Nut-KO_over_Input-Nut-WT.pdf
        # Same with subsampling:
        out/r/frip_from_macs2_output/macs2/callpeak_broad_format-BAMPE_gsize-mm_SPMR/samtools/view/subsampling_bam/seed7_30000000pf/picard/MarkDuplicates/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run140_run141/H4K5ac-Nut-WT_over_Input-Nut-WT_VS_Input-Nut-KO_over_Input-Nut-WT.pdf

    """
    input:
        xls1="out/{filler}/{id1}_peaks.xls",
        xls2="out/{filler}/{id2}_peaks.xls",
        code="src/r/script/frip_from_macs2_output.R"
    output:
        pdf="out/r/frip_from_macs2_output/{filler}/{id1}_VS_{id2}.pdf"
    conda:
        "../envs/r.yaml"
    shell:
        "Rscript {input.code} -i {input.xls1} -j {input.xls2} -o {output.pdf}"
    #run:
    #    R("""
    #    library(beanplot)

    #    xls1 <- read.table(
    #        #file='out/macs2/callpeak_broad_BAMPE_mm_SPMRTrue/bam/mm9/merge_run140_run141/H4K5ac-Nut-WT_over_Input-Nut-WT_peaks.xls',
    #        file='{input.xls1}',
    #        comment.char = "#",
    #        header = T)

    #    xls2 <- read.table(
    #        #file='out/macs2/callpeak_broad_BAMPE_mm_SPMRTrue/bam/mm9/merge_run140_run141/H4K5ac-Nut-KO_over_Input-Nut-KO_peaks.xls',
    #        file='{input.xls2}',
    #        comment.char = "#",
    #        header = T)

    #    # Having trouble sometimes with default bw settings:
    #    # http://r.789695.n4.nabble.com/beanplot-Error-sample-is-too-sparse-to-find-TD-td4291739.html

    #    pdf('{output.pdf}')
    #    beanplot(
    #        xls1$pileup,
    #        xls2$pileup,
    #        xls1$length,
    #        xls2$length,
    #        bw="nrd0", 
    #        side='both',
    #        names=c('pileup','length'),
    #        ylab='MACS2 metrics (unit in class labels)',
    #        xlab='WT (half left) versus KO (half right)',
    #        main=paste0('n peaks WT: ', dim(xls1)[1], '; KO: ', dim(xls2)[1], '.')
    #        )
    #    dev.off()
    #    """)


