rule bedtools_reldist:
    """
    Created:
        2017-12-08 15:43:09
    Test:
        features_a:
            out/awk/extract_xls_coordinates_to_bed3/macs2/callpeak_format-BAMPE_gsize-mm_SPMR/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/sickle/pe_q-20_t-sanger/ln/pe_remove_mate_prefix/gunzip/merge_lanes_nextseq500_pe/ln/alias/fastq/run200/P300-WT_over_input-WT_peaks.bed
        features_b:
            out/awk/extract_xls_coordinates_to_bed3/macs2/callpeak_format-BAMPE_gsize-mm_SPMR/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/sickle/pe_q-20_t-sanger/ln/pe_remove_mate_prefix/gunzip/merge_lanes_nextseq500_pe/ln/alias/fastq/run200/Nut-WT_over_input-WT_peaks.bed
        output:
            out/bedtools/reldist/b-awk/extract_xls_coordinates_to_bed3/macs2/callpeak_format-BAMPE_gsize-mm_SPMR/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/sickle/pe_q-20_t-sanger/ln/pe_remove_mate_prefix/gunzip/merge_lanes_nextseq500_pe/ln/alias/fastq/run200/P300-WT_over_input-WT_peaks_bed/a-awk/extract_xls_coordinates_to_bed3/macs2/callpeak_format-BAMPE_gsize-mm_SPMR/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/sickle/pe_q-20_t-sanger/ln/pe_remove_mate_prefix/gunzip/merge_lanes_nextseq500_pe/ln/alias/fastq/run200/P300-WT_over_input-WT_peaks_bed.tsv

        next step
    """
    input:
        features_a="out/{filler1}.{ext1}",
        features_b="out/{filler2}.{ext2}"
    output:
        tsv="out/bedtools/reldist/b-{filler2}_{ext2}/a-{filler1}_{ext1}.tsv"
    wildcard_constraints:
        ext1="bam|bed|bedgraph|gff|vcf",
        ext2="bam|bed|bedgraph|gff|vcf"
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        bedtools\
            reldist\
            -a {input.features_a}\
            -b {input.features_b}\
            > {output.tsv}
        """

