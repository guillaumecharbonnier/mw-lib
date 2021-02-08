rule rseg_deadzones_kmer_prefix_suffix:
    """
    Created:
        2017-10-31 14:26:25
    Aim:
        The deadzone program in RSEG software package is used to compute unmappable regions given genome
        assembly and read length.  You need first to download the genome sequence of the genome in fasta format
        from UCSC Genome Browser Download. Suppose the fasta files containing the sequence for mouse mm9 is
        located at mm9/. You can compute unmappable regions for 32bp reads by running the following command.
    Doc:
        -o, -output   Name of output file (default: stdout) 
        -k, -kmer     Width of k-mers 
        -p, -prefix   prefix length (default 5) 
        -s, -suffix   suffix of FASTA files (assumes -c indicates dir) 
        -v, -verbose  print more run information
    Note:
        out/gunzip/to-stdout/wget/ftp_ensembl/pub/release-87/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.chromosome.1.fa
        out/cat/assembly_ensembl/GRCm38.fa
        #filename={input.fasta}
        #filenameNoSuffix="${{filename%.*}}"
    Test:
        out/rseg/deadzones_kmer-32_prefix-5_suffix-fa/cat/assembly_ensembl/GRCm38.bed
        out/rseg/deadzones_kmer-32_prefix-5_suffix-fa/cat/assembly_ensembl/GRCh38.bed
    """
    input:
        fasta="out/{filler}.{suffix}"
    output:
        bed="out/rseg/deadzones_kmer-{kmer}_prefix-{prefix}_suffix-{suffix}/{filler}.bed"
    wildcard_constraints:
        kmer="[0-9]+",
        prefix="[0-9]+",
        suffix="fa|fasta",
    conda:
        "../envs/rseg.yaml"
    shell:
        """
        deadzones \
            -output {output.bed} \
            -kmer {wildcards.kmer} \
            -prefix {wildcards.prefix} \
            -suffix {wildcards.suffix} \
            {input.fasta}
        """

rule rseg_diff_chrom_deadzones_mode:
    """
    Created:
        2017-10-31 15:17:31
    Aim:
        rseg-diff can be used in two ways: first, it is used to find histone domains by using both a test sample and a
        control sample.  Second, it is used to find domains with different signals either between two histone marks
        in the same cell type or between two cell types with the same histone modifcation.
    Warning:
        rseg requires the input read files are sorted by chromosome name, starting position, ending position and strand, which can be done with standard UNIX sort tool as following:
        $ export LC_ALL=C
        $ sort -k1,1 -k2,2n -k3,3n -k6,6 -o sorted.bed input.bed
        Note that we need to set the locale of the shell environment to the C programming language locale.
    Doc:
        Usage: rseg-diff [OPTIONS] <mapped-read-locations-A> <mapped-read-locations-B>
        -o, -out               domain output file 
            -score             Posterior scores file 
            -readcount         readcounts file 
            -boundary          domain boundary file 
            -boundary-score    boundary transition scores file 
        -c, -chrom             file with chromosome sizes (BED format) 
        -d, -deadzones         file of deadzones (BED format) 
        -B, -bam               Input reads file is BAM format 
            -param-in          Input parameters file 
            -param-out         Output parameters file 
        -m, -mode              running mode 2:test-control; 3: test-test 
        -i, -maxitr            maximum iterations for training 
        -b, -bin-size          bin size (default: based on data) 
            -bin-step          minimum bin size (default: 50) 
            -duplicates        keep duplicate reads 
            -fragment_length   Extend reads to fragment length (default not to 
                               extend) 
            -Waterman          use Waterman's method for bin size 
            -Hideaki           use Hideaki's method for bin size 
            -Hideaki-emp       use Hideaki's empirical method (default) 
            -smooth            Indicate whether the rate curve is assumed smooth 
            -max-dead          max deadzone proportion for retained bins 
        -s, -domain-size       expected domain size (default: 20000) 
        -S, -desert            desert size (default: 20000) 
        -F, -fg                foreground emission distribution 
        -B, -bg                background emission distribution 
            -training-size     Max number of data points for training (default: all) 
        -P, -posterior         use posterior decoding (default: Viterbi) 
            -posterior-cutoff  Posterior threshold for signigicant bins 
            -undefined         min size of unmappable region 
            -cutoff            cutoff in cdf for identified domains 

    Test:
        Previous test using bam mode:
            out/rseg_diff_chrom-GRCm38_deadzones-GRCm38_bam_mode-3/samtools/sort/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_GRCm38/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run140_run141/H4K5ac-Nut-WT_VS_H4K5ac-Nut-KO.bed 
            out/rseg_diff_chrom-mm9_deadzones-mm9_bam_mode-3/samtools/sort/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm9/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run140_run141/H4K5ac-Nut-WT_VS_H4K5ac-Nut-KO.bed 

        out/rseg_diff_chrom-GRCm38_deadzones-GRCm38_mode-3/sort/coordinates_rseg_bed/bedtools/bamtobed/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_GRCm38/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run140_run141/H4K5ac-Nut-WT_VS_H4K5ac-Nut-KO.bed
        out/rseg_diff_chrom-mm9_deadzones-mm9_mode-3/sort/coordinates_rseg_bed/bedtools/bamtobed/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm9/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run140_run141/H4K5ac-Nut-WT_VS_H4K5ac-Nut-KO.bed
        Only chr9 to test on small dataset:
            out/rseg_diff_chrom-mm9_deadzones-mm9_mode-3/awk/extract_chr9/sort/coordinates_rseg_bed/bedtools/bamtobed/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm9/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run140_run141/H4K5ac-Nut-WT_VS_H4K5ac-Nut-KO.bed

        Use pairedend reads bampe:
            out/rseg_diff_chrom-mm9_deadzones-mm9_mode-3/sort/coordinates_rseg_bed/awk/convert_bedpe_to_bed6_insert_size/bedtools/bamtobed_bedpe/samtools/sort_n/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm9/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run140_run141/H4K5ac-Nut-WT_VS_H4K5ac-Nut-KO.bed
            out/rseg_diff_chrom-mm9_deadzones-mm9_mode-3/sort/coordinates_rseg_bed/awk/convert_bedpe_to_bed6_insert_size/bedtools/bamtobed_bedpe/samtools/sort_n/samtools/merge/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_mm9/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run140_run141/H4K5ac-Nut-WT_VS_H4K5ac-Nut-KO.bed
        Test for nucleosomes and small structures:
            First test:
                out/rseg_diff_chrom-mm9_deadzones-mm9_mode-3/sort/coordinates_rseg_bed/awk/convert_bedpe_to_bed6_insert_size/bedtools/bamtobed_bedpe/samtools/sort_n/samtools/merge/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_mm9/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band1_VS_MNase_Spm_WT_band2.bed
            Better test with all runs and size selection:
                input:
                    spermatozoa_mnase_nuc="out/awk/extract_main_chr/awk/convert_bedpe_to_bed6_insert_size/bedtools/bamtobed_bedpe/samtools/sort_n/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-130_lmax-170/samtools/merge_three_samples/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167_run184_run187/s1-MNase_Spm_WT_band1_s2-MNase_Spm_WT_band3_s3-MNase_Spm_WT_band6.bed"
                    spermatozoa_mnase_ss="out/awk/extract_main_chr/awk/convert_bedpe_to_bed6_insert_size/bedtools/bamtobed_bedpe/samtools/sort_n/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-30_lmax-100/samtools/merge_three_samples/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167_run184/s1-MNase_Spm_WT_band2_s2-MNase_Spm_WT_band4_s3-MNase_Spm_WT_band5.bed"
                output:
                    out/rseg/diff_chrom-mm9_deadzones-mm9_mode-3/sort/coordinates_rseg_bed/ln/alias/mm9/spermatozoa_mnase_nuc_VS_spermatozoa_mnase_ss.bed
                    out/rseg/diff_chrom-mm9_deadzones-mm9_mode-3/sort/coordinates_rseg_bed/ln/alias/mm9/spermatozoa_mnase_nuc_VS_round_mnase_nuc.bed
    Test:
        out/rseg/diff_chrom-GRCh38-chrom_deadzones-GRCh38-deadzones_mode-3/sort/coordinates-rseg-bed/bedtools/bamtobed/ln/alias/experiments/hg38_H3K27ac_thymus/CD34_VS_EC.bed

    Note:
        -chrom can't be read if put before deadzones? not sure but strange. Actually it seems that the parser does not like multiple spaces between arguments and the command in one line works correctly.
    """
    input:
        chrom     = lambda wildcards: eval(mwconf['ids'][wildcards.chrom_id]),
        deadzones = lambda wildcards: eval(mwconf['ids'][wildcards.bed_id]),
        bed1 = "out/{filler}/{sample1}.bed",
        bed2 = "out/{filler}/{sample2}.bed"
    output:
        #bed            = "out/rseg/diff_chrom-{chrom_id}_deadzones-{bed_id}_mode-{mode}/{filler}/{sample1}_VS_{sample2}.bed",
        #score          = "out/rseg/diff_chrom-{chrom_id}_deadzones-{bed_id}_mode-{mode}/{filler}/{sample1}_VS_{sample2}.score.txt",
        #readcount      = "out/rseg/diff_chrom-{chrom_id}_deadzones-{bed_id}_mode-{mode}/{filler}/{sample1}_VS_{sample2}.readcount.txt",
        #boundary       = "out/rseg/diff_chrom-{chrom_id}_deadzones-{bed_id}_mode-{mode}/{filler}/{sample1}_VS_{sample2}.boundary.txt",
        #boundary_score = "out/rseg/diff_chrom-{chrom_id}_deadzones-{bed_id}_mode-{mode}/{filler}/{sample1}_VS_{sample2}.boundary_score.txt",
        bed            = "out/rseg/diff_{chrom_id}_deadzones-{bed_id}_mode-{mode}/{filler}/{sample1}_VS_{sample2}.bed",
        score          = "out/rseg/diff_{chrom_id}_deadzones-{bed_id}_mode-{mode}/{filler}/{sample1}_VS_{sample2}.score.txt",
        readcount      = "out/rseg/diff_{chrom_id}_deadzones-{bed_id}_mode-{mode}/{filler}/{sample1}_VS_{sample2}.readcount.txt",
        boundary       = "out/rseg/diff_{chrom_id}_deadzones-{bed_id}_mode-{mode}/{filler}/{sample1}_VS_{sample2}.boundary.txt",
        boundary_score = "out/rseg/diff_{chrom_id}_deadzones-{bed_id}_mode-{mode}/{filler}/{sample1}_VS_{sample2}.boundary_score.txt",
    wildcard_constraints:
        mode = "2|3",
        chrom_id = "[a-zA-Z0-9-]+",
        bed_id = "[a-zA-Z0-9-]+"
    conda:
        "../envs/rseg.yaml"
    shell:
        """
        rseg-diff\
            -out {output.bed}\
            -score {output.score}\
            -readcount {output.readcount}\
            -boundary {output.boundary}\
            -boundary-score {output.boundary_score}\
            -deadzones {input.deadzones}\
            -chrom {input.chrom}\
            -mode {wildcards.mode}\
            -v\
            {input.bed1} {input.bed2}
        """

