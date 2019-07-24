"""
Created:
    2017-03-08 10:45:14
Aim:
Doc:
    http://crossmap.sourceforge.net/#usage
    Program:
        CrossMap (v0.2.2)

    Description: 
        CrossMap is a program for convenient conversion of genome coordinates and
        genome annotation files between assemblies (eg. lift from human hg18 to hg19 or
        vice versa).It supports file in BAM, SAM, BED, Wiggle, BigWig, GFF, GTF and VCF
        format.

    Usage:
        CrossMap.py <command> [options]

        bam   convert alignment file in BAM or SAM format.
        bed   convert genome cooridnate or annotation file in BED or BED-like format.
        bigwig        convert genome coordinate file in BigWig format.
        gff   convert genome cooridnate or annotation file in GFF or GTF format.
        vcf   convert genome coordinate file in VCF format.
        wig   convert genome coordinate file in Wiggle, or bedGraph format.
"""
rule crossmap_bed_sam_gff_gtf_vcf:
    """
    Created:
        2017-03-08 11:06:16
    Aim:
        I need to convert BRDT peaks from mm9 to mm10 to be able to use gtftk coverage which requires a complete gtf with 'gene' and 'transcript'.
    Warning:
        Crossmap may not work (not tested) for the listed format in the rule name (see rules for specific formats below).
    Note:
        Why bedtools merge may be needed after crossmapping:
        $ wc -l inp/bed/atac/cd34_sort.bed out/crossmap/hg38_to_hg19/inp/bed/atac/cd34_sort.bed out/bedtools/merge_-d_2/sort/coordinates_bed/crossmap/hg38_to_hg19/inp/bed/atac/cd34_sort.bed out/bedtools/merge_-d_5/sort/coordinates_bed/crossmap/hg38_to_hg19/inp/bed/atac/cd34_sort.bed out/bedtools/merge_-d_10/sort/coordinates_bed/crossmap/hg38_to_hg19/inp/bed/atac/cd34_sort.bed out/bedtools/merge_-d_20/sort/coordinates_bed/crossmap/hg38_to_hg19/inp/bed/atac/cd34_sort.bed out/bedtools/merge_-d_15/sort/coordinates_bed/crossmap/hg38_to_hg19/inp/bed/atac/cd34_sort.bed out/bedtools/merge_-d_12/sort/coordinates_bed/crossmap/hg38_to_hg19/inp/bed/atac/cd34_sort.bed out/bedtools/merge_-d_13/sort/coordinates_bed/crossmap/hg38_to_hg19/inp/bed/atac/cd34_sort.bed out/bedtools/merge_-d_14/sort/coordinates_bed/crossmap/hg38_to_hg19/inp/bed/atac/cd34_sort.bed
        2522 inp/bed/atac/cd34_sort.bed
        2541 out/crossmap/hg38_to_hg19/inp/bed/atac/cd34_sort.bed
        2526 out/bedtools/merge_-d_2/sort/coordinates_bed/crossmap/hg38_to_hg19/inp/bed/atac/cd34_sort.bed
        2525 out/bedtools/merge_-d_5/sort/coordinates_bed/crossmap/hg38_to_hg19/inp/bed/atac/cd34_sort.bed
        2524 out/bedtools/merge_-d_10/sort/coordinates_bed/crossmap/hg38_to_hg19/inp/bed/atac/cd34_sort.bed
        2519 out/bedtools/merge_-d_20/sort/coordinates_bed/crossmap/hg38_to_hg19/inp/bed/atac/cd34_sort.bed
        2520 out/bedtools/merge_-d_15/sort/coordinates_bed/crossmap/hg38_to_hg19/inp/bed/atac/cd34_sort.bed
        2523 out/bedtools/merge_-d_12/sort/coordinates_bed/crossmap/hg38_to_hg19/inp/bed/atac/cd34_sort.bed
        2523 out/bedtools/merge_-d_13/sort/coordinates_bed/crossmap/hg38_to_hg19/inp/bed/atac/cd34_sort.bed
        2523 out/bedtools/merge_-d_14/sort/coordinates_bed/crossmap/hg38_to_hg19/inp/bed/atac/cd34_sort.bed

    Todo:
        I should add exception for bigwig when ext is 'bw'.
    Note:
        Because output file is not always the same type as input file. I split crossmap rule into the general rule for formats where they are the same. And add below another rule for exceptions.
        NOT TESTED FOR ALL FORMATS.
        Output file
        
        Format of Output files depends on the input format
        Input_format    Output_format
        BED     BED (Genome coordinates will be updated to the target assembly)
        BAM     BAM (Genome coordinates, header section, all SAM flags, insert size will be updated accordingly)
        SAM     SAM (Genome coordinates, header section, all SAM flags, insert size will be updated accordingly)
        Wiggle  bedGraph (if wigToBigWig executable does not exist)
        Wiggle  BigWig (if wigToBigWig executable exists)
        BigWig  bedGraph (if wigToBigWig executable does not exist)
        BigWig  BigWig (if wigToBigWig executable exists)
        GFF     GFF (Genome coordinates will be updated to the target assembly)
        GTF     GTF (Genome coordinates will be updated to the target assembly)
        VCF     VCF (Genome coordinates and reference alleles will be updated to the target assembly)
   Test:
        out/crossmap/bed_hg19_to_hg38/inp/bed/salva/Regions_capture_SE_Starr-seq_Alex_HG19_Merged_cellLine_specific_SE.bed
        out/crossmap/bed_mm9_to_mm10/bedtools/intersect/tss_with_brdt_and_h4k5ac.bed
        out/crossmap/bed_hg38_to_hg19/inp/2018_02_22_cluster10_debug_thymus_to_great.bed
        out/crossmap/bed_hg38_to_hg19/awk/extract_main_chr/bedtools/bamtobed/inp/bam/GRCh38/Blueprint/TH89_EC_S00FJ6H1.ERX1302276.H3K27ac.dedup.bwa.GRCh38.20160225.bed
        out/ucsc/wigToBigWig_hg38-main-chr_clip/crossmap/hg19_to_hg38/gunzip/to-stdout/wget/https/www.genboree.org/EdaccData/Current-Release/experiment-sample/Bisulfite-Seq/Mobilized_CD34_Primary_Cells/BI.Mobilized_CD34_Primary_Cells.Bisulfite-Seq.RO_01549.bw
        out/crossmap/hg19_to_hg38/gunzip/to-stdout/wget/https/www.genboree.org/EdaccData/Current-Release/experiment-sample/Bisulfite-Seq/Mobilized_CD34_Primary_Cells/BI.Mobilized_CD34_Primary_Cells.Bisulfite-Seq.RO_01549.wig
        out/crossmap/bed_hg38_to_hg19/r/dynamic_enchancers_in_thymopoiesis/dClust/rowFeature-no_rmsk_mxy__no_donor_effect__distal/sortingMethod-kmeans_centers-8_nstart-100_itermax-200000_algorithm-Lloyd/sortingFeature-cpg_meth_call/subSortingMethod-kmeans_centers-8_nstart-100_itermax-200000_algorithm-Lloyd/subSortingFeature-H3K27ac_peaks.ATAC_peaks/cluster-1.bed
        out/crossmap/chain-mm10-to-mm9/cut/_-f1-6/macs2/callpeak_--broad/samtools/index/samtools/sort/samtools/view_sam_to_bam_-q_30/bowtie2/se_mm10/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR3126243_over_SRR3126242_peaks.bed
    """
    input:
        chain = lambda wildcards: config['ids'][wildcards.chain_id],
        coord = "out/{filler}.{ext}"
    output:
        coord = "out/crossmap/{chain_id}/{filler}.{ext}"
    wildcard_constraints:
        ext = "bed|bigwig|gff|vcf|wig"
    conda:
        "../envs/crossmap.yaml"
    shell:
        "CrossMap.py {wildcards.ext} {input.chain} {input.coord} {output.coord}"

rule crossmap_wig:
    """
    Created:
        2017-03-08 11:06:16
    Aim:
        Specific rule for wig files.
    Test:
        out/crossmap/wig_hg19_to_hg38/gunzip/to-stdout/wget/https/www.genboree.org/EdaccData/Current-Release/experiment-sample/Bisulfite-Seq/Mobilized_CD34_Primary_Cells/BI.Mobilized_CD34_Primary_Cells.Bisulfite-Seq.RO_01549.bw
    """
    input:
        chain = lambda wildcards: config['ids'][wildcards.chain_id],
        coord = "out/{filler}.wig"
    output:
        coord = "out/crossmap/wig_{chain_id}/{filler}.bw"
    conda:
        "../envs/crossmap.yaml"
    shell:
        """
        CrossMap.py wig {input.chain} {input.coord} {output.coord}
        # Because crossmap expects to be given the prefix and not the output bw filename:
        mv {output.coord}.bw {output.coord}
        """

rule crossmap_bam:
    """
    Created:
        2018-08-20 15:14:41
    Aim:
        I need to convert Solid Bam from hg19 to hg38 for the thymus project.
    Test:
    """
    input:
        chain = lambda wildcards: config['ids'][wildcards.chain_id],
        bam="out/{filler}.bam"
    output:
        bam="out/crossmap/{chain_id}/{filler}.bam",
        bai="out/crossmap/{chain_id}/{filler}.bam.bai",
    conda:
        "../envs/crossmap.yaml"
    shell:
        """
        CrossMap.py bam {input.chain} {input.bam} {output.bam}
        mv {output.bam}.sorted.bam {output.bam} 
        mv  {output.bam}.sorted.bam.bai {output.bai}
        # crossmap let the unsorted bam behind:
        rm -f {output.bam}.bam
        """


