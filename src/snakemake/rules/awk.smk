rule awk_extract_rows_matching_col_content:
    """
    out/awk/extract_rows_matching_col13_content/deepTools/plotHeatmap_--kmeans_40_--boxAroundHeatmaps_no/deepTools/computeMatrix_reference-point_--referencePoint_center_-bs_5_-b_150_-a_150_bed-mm10-nuc-prcs-maxfuzz-40-minsmt-30-only-high-qual-fragments_bw-mnase-prcs-naked-nuc-ss/cluster_1.bed
    """
    input:
        "out/{filler}.bed"
    output:
        "out/awk/extract_rows_matching_col{col}_content/{filler}/{content}.bed"
    wildcard_constraints:
        col="[0-9]+",
        content="[\w_-]+"
    shell:
        """
        awk '${wildcards.col} == "{wildcards.content}"' {input} > {output}
        """

rule awk_trim_sam_query_name:
    """
    Created:
        2019-11-29 11:14:07
    Aim:
        SAM format limit query name to 254 characters.
        Sometimes this query can overcome this limit, 
        especially using Metaworkflow paradigm.
        This rule trim the query to make sam file
        in compliance with the standard.
    """
    input:
        "out/{filler}.sam"
    output:
        "out/awk/trim_sam_query_name/{filler}.sam"
    shell:
        """
        awk 'BEGIN{{OFS=FS="\\t"}} {{$1=substr($1,1,250)}}1' {input} > {output}
        """

rule awk_extract_cluster:
    """
    out/awk/extract_cluster/deepTools/plotHeatmap_--kmeans_40_--boxAroundHeatmaps_no/deepTools/computeMatrix_reference-point_--referencePoint_center_-bs_5_-b_150_-a_150_bed-mm10-nuc-prcs-maxfuzz-40-minsmt-30-only-high-qual-fragments_bw-mnase-prcs-naked-nuc-ss/1.bed
    """
    input:
        "out/{filler}.bed"
    output:
        "out/awk/extract_cluster/{filler}/{cluster}.bed"
    wildcard_constraints:
        cluster="[0-9]+"
    shell:
        """
        awk '$13 == "cluster_{wildcards.cluster}"' {input} > {output}
        """

rule awk_collapse_Carrillo2017_11_states_to_5:
    """
    Aim:
        Take a segmentation file from ChromHMM and collapse the 11 states into 2 states
    Note:
        Segmentation model correspond to:
        1   T
        2   T
        3   H
        4   H
        5   H
        6   H
        7   R
        8   E
        9   E
        10  A
        11  A
    Test:
        out/awk/collapse_Carrillo2017_11_states_to_5/ChromHMM/MakeSegments_blueprint_tall_samples_into_chromdet_original_paper/T11C_11_segments.bed
    """
    input:
        "out/{filler}.bed"
    output:
        "out/awk/collapse_Carrillo2017_11_states_to_5/{filler}.bed"
    shell:
        """
        awk '
            BEGIN{{FS=OFS="\\t"}}
            {{
            if ($4 == "E1" || $4 == "E2") {{ $4 = "T" }}
            if ($4 == "E3" || $4 == "E4" || $4 == "E5" || $4 == "E6") {{ $4 = "H" }}
            if ($4 == "E7") {{ $4 = "R" }}
            if ($4 == "E8" || $4 == "E9") {{ $4 = "E" }}
            if ($4 == "E10" || $4 == "E11") {{ $4 = "A" }}
            print $0 }}
        ' {input} > {output}
        """

rule awk_merge_book_ended_features_in_bed_if_same_name:
    """
    From:
        https://bioinformatics.stackexchange.com/questions/3523/bed-file-merge-book-end-features-only-if-same-name-in-column-4
    Test:
        out/awk/merge_book_ended_features_in_bed_if_same_name/awk/collapse_Carrillo2017_11_states_to_5/ChromHMM/MakeSegments_blueprint_tall_samples_into_chromdet_original_paper/T11C_11_segments.bed
    """
    input:
        "out/{filler}"
    output:
        "out/awk/merge_book_ended_features_in_bed_if_same_name/{filler}"
    shell:
        """
        awk 'BEGIN{{OFS="\t"}}
        {{
            if($1==lchrom && $4==lname && $2 == lend) {{lend = $3}}
            else{{
                if(lchrom) {{print lchrom, lstart, lend, lname;}}; lchrom=$1; lstart=$2; lend=$3; lname=$4}}
        }}END{{print lchrom, lstart, lend, lname}}' {input} > {output}
        """

rule awk_tfbsConsSites_to_gtf:
    """
    """
    input:
        "out/gunzip/_-c/rsync/_-aP/hgdownload.cse.ucsc.edu/goldenPath/hg19/database/tfbsConsSites.txt"
    output:
        "out/awk/tfbsConsSites_to_gtf.gtf"
    log:
        "out/awk/tfbsConsSites_to_gtf.log"
    benchmark:
        "out/awk/tfbsConsSites_to_gtf.benchmark.tsv"
    shell:
        """
        awk -v q='"' 'BEGIN{{FS=OFS="\\t"}} {{gsub("V.", "", $5); gsub("_.+$", "", $5); print $2, "tfbsConsSites", "gene", $3, $4, $8, $7, ".", "gene_id " q NR q "; tf_id " q $5 q}}' {input} > {output} 2> {log}
        """

rule awk_fix_bed9_thick_cols:
    """
    Created:
        2019-10-13 23:39:19
    Aim:
        When using crossmap, thick start and end coordinates are strangely converted.
        For segmentation tracks, I only need to have the same coordinates as start and end so I replicate these columns.
    """
    input:
        bed="out/{filler}.bed"
    output:
        bed="out/awk/fix_bed9_thick_cols/{filler}.bed"
    shell:
        """
        awk 'BEGIN{{FS=OFS="\\t"}} {{print $1,$2,$3, $4, $5, $6, $2, $3, $9}}' {input.bed} > {output.bed}
        """

rule awk_rename_bed_name:
    """
    Created:
        2017-01-24 18h14
    Test:
        out/awk/rename_bed_name/sort/unique_coordinates_bed/input/annotation/feature/debugged_ranked_florent/mm9_5000_features_dwnreg.bed
    """
    input:
        bed="out/{filler}.bed"
    output:
        bed="out/awk/rename_bed_name/{filler}.bed"
    shell:
        """
        awk 'BEGIN{{FS=OFS="\\t"}} {{print $1,$2,$3, NR, $5, $6}}' {input.bed} > {output.bed}
        """

rule awk_remove_bed_name:
    """
    Created:
        2019-01-14 15:42:57
    Aim:
        Replace bed name (4th column) by '.'.
        This is useful when looking at ChromHMM segmentation files in IGV because we do not care about state number and we just want to see the color.
    Test:
        out/awk/remove_bed_name/crossmap/hg19_to_hg38/ChromHMM/MakeBrowserFiles/ChromHMM/MakeSegments_our_merged_thymic_samples_into_chromdet_original_paper/CD34_11_segments_dense.bed
    """
    input:
        bed="out/{filler}.bed"
    output:
        bed="out/awk/remove_bed_name/{filler}.bed"
    shell:
        """
        awk 'BEGIN{{FS=OFS="\\t"}} {{$4="."; print}}' {input.bed} > {output.bed}
        """


rule awk_fill_bed3_to_bed6:
    """
    Created:
        2017-03-15 17:35:12
    Aim:
        I have bed3 and I need to fill name, strand and score for GFOLD to deal with it as annotation reference.
    Test:
        out/awk/fill_bed3_to_bed6/bedtools/merge/sort/coordinates_bed/cat/mm10_nut_wt_ko_h4k5ac_peaks.bed
    """
    input:
        bed="out/{filler}.bed"
    output:
        bed="out/awk/fill_bed3_to_bed6/{filler}.bed"
    shell:
        """
        awk 'BEGIN{{FS=OFS="\\t"}} {{print $1,$2,$3, $1"-"$2"-"$3, $3-$2, "."}}' {input.bed} > {output.bed}
        """

rule awk_extract_main_chr:
    """
    Created:
        2017-03-08 14:05:37
    Test:
        "out/awk/extract_main_chr/gunzip/rsync/ucsc/goldenPath/mm10/database/chromInfo.txt"
    """
    input:
        txt="out/{filler}"
    output:
        txt="out/awk/extract_main_chr/{filler}"
    shell:
        """
        awk '/^chr[0-9XYM]+\t/' {input.txt} > {output.txt}
        """

rule awk_extract_autosomes:
    """
    Created:
        2018-04-30 11:38:15
    Test:
        "out/awk/extract_autosomes/gunzip/rsync/ucsc/goldenPath/mm10/database/chromInfo.txt"
    """
    input:
        txt="out/{filler}"
    output:
        txt="out/awk/extract_autosomes/{filler}"
    shell:
        """
        awk '/^chr[0-9]+\t/' {input.txt} > {output.txt}
        """

rule awk_extract_chr:
    """
    Created:
        2018-03-07 15:51:06
    Test:
        "out/awk/extract_chr1/gunzip/rsync/ucsc/goldenPath/mm10/database/chromInfo.txt"
        out/awk/extract_chr1/sed/add_chr/gunzip/to-stdout/wget/ftp/ftp.ensembl.org/pub/release-91/gff3/mus_musculus/Mus_musculus.GRCm38.91.gff3
    """
    input:
        txt="out/{filler}"
    output:
        txt="out/awk/extract_chr{chr}/{filler}"
    wildcard_constraints:
        chr="[0-9]+|X|Y|M"
    shell:
        """
        awk '/^chr{wildcards.chr}+\t/' {input.txt} > {output.txt}
        """


rule awk_extend_reads:
    """
    Modified:
        2017-04-11 14:36:26 - Integrated to this workflow
    Note:
        Taken from Aurelien's workflow.
        Added double escape for tabulation in awk.
        Added check condition if chromosomes are like "chr3" or only "3".
    """
    input:
        bed="out/{filler}.bed"
    output:
        bed="out/awk/extend_reads_{extReads}/{filler}.bed"
    wildcard_constraints:
        extReads="[0-9]+"
    shell:
        """
        size_fragment={wildcards.extReads}
        awk -v FRAG_SIZE=$size_fragment 'BEGIN{{ OFS="\\t" }}{{
            if ($1 ~ "^chr*"){{
                if ($6 == "+"){{ $3 = $2 + FRAG_SIZE ; print $0 }}
                else if ($6 == "-"){{ $2 = $3 - FRAG_SIZE ; if ($2 > 0) {{ print $0 }} }}
                }}
            else {{
                if ($6 == "+"){{ $3 = $2 + FRAG_SIZE ; print "chr"$0 }}
                else if ($6 == "-"){{ $2 = $3 - FRAG_SIZE ; if ($2 > 0) {{ print "chr"$0 }} }}
                }}
            }}' {input.bed} > {output.bed}
        """

rule awk_keep_first_mate_for_pe_bedtools_bamtobed:
    """
    Created:
        2018-01-30 12:35:22
    Aim:
        Keep only the first mate of paired-end reads when they are converted from bam to bed by bedtools bamtobed.
        It works based on the pattern naming of both mates:
        chr7   118970079   118970129   TUPAC_0001:3:1:0:1452#0/1   37   -
        chr7   118965072   118965122   TUPAC_0001:3:1:0:1452#0/2   37   +
        The first use for this rule is to merge samples which are se and pe for CapStarrSeq.
    Test:

    """
    input:
        bed="out/{filler}.bed"
    output:
        bed="out/awk/keep_first_mate_for_pe_bedtools_bamtobed/{filler}.bed"
    shell:
        """
        awk 'BEGIN{{ OFS="\\t" }}{{
            if ($4 !~ "/2$"){{print $0 }}
            }}' {input.bed} > {output.bed}
        """


rule awk_chromInfo_to_bed3:
    """
    Created:
        2017-04-12 14:39:46
    Modified:
        2018-08-25 18:22:17
    Aim:
        I need to produce a bedfile for bedtools_complementary_features from a chrominfo.
    """
    input:
        "out/{filler}.txt"
    output:
        "out/awk/chromInfo_to_bed3/{filler}.bed"
    shell:
        """
        awk 'BEGIN{{FS=OFS="\\t" }} {{print $1,0,$2}}' {input} > {output}
        """

rule awk_faidx_to_bed3:
    """
    Created:
        2017-11-02 11:20:07
    Aim:
        I need to produce a chromInfo-like bedfile for RSEG using Ensembl GRCm38.
    Test:
        out/awk/faidx_to_bed3/cat/assembly_ensembl/GRCm38.bed
    """
    input:
        "out/{filler}.fa.fai"
    output:
        "out/awk/faidx_to_bed3/{filler}.bed"
    shell:
        """
        awk 'BEGIN{{FS=OFS="\\t" }} {{print $1,1,$2}}' {input} > {output}
        """

rule awk_extract_two_first_columns:
    """
    Created:
        2017-04-12 15:22:50
        extract_two_first_columns
    Aim:
        Bedtools complement does not like the fact that UCSC chromInfo have a third useless column.
    """
    input:
        "out/{filler}"
    output:
        "out/awk/extract_two_first_columns/{filler}"
    shell:
        """
        awk 'BEGIN{{FS=OFS="\\t" }} {{print $1,$2}}' {input} > {output}
        """

rule awk_extract_danpos_pos_to_bed_maxfuz_minsmt:
    """
    Created:
        2017-02-06 10:16:39
    Aim:
        Extract positions with low fuzziness for small structures in order to keep only those correctly called by Danpos and use them to look if they are enriched in some motif in the center. I could also intersect these selected positions with the nucleosomes positions with low fuzziness so I look only at regions that were containing stable nucleosomes in round spermatid.
    Test:
        out/awk/extract_danpos_pos_to_bed_maxfuz-50_minsmt-100/danpos/dtriple_v4/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run113_run125_run126/MNS-SC-WT/danpos.smooth.positions.bed
    """
    input:
        xls="out/{filler}.xls"
    output:
        bed="out/awk/extract_danpos_pos_to_bed_maxfuz-{maxfuz}_minsmt-{minsmt}/{filler}.bed"
    wildcard_constraints:
        maxfuz="[0-9]+",
        minsmt="[0-9]+"
    shell:
        """
        awk 'BEGIN {{OFS="\\t"}}; ($5 > {wildcards.minsmt}) && ($6 < {wildcards.maxfuz}) {{print $1, $2, $3;}}' {input.xls} > {output.bed}
        """

rule awk_shrink_bed:
    """
    Created:
        2016-12-12 13h35
    Aim:
        This rule does the opposite of Bedtools slop and reduces the range of coordinates.
        It is used to resize danpos positions to small structure size, shrinking from 140bp to 70 with cutbp=35.
    Test:
        input:"out/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126/PSK-SC-WT/danpos.smooth.positions.bed"
        output:"out/awk/shrink_bed/35bp/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126/PSK-SC-WT/danpos.smooth.positions.bed"
    """
    input:
        bed="out/{filler}.bed"
    output:
        bed="out/awk/shrink_bed/{cutbp}bp/{filler}.bed"
    wildcard_constraints:
        cutbp="[0-9]+"
    shell:
        """
        awk 'BEGIN {{OFS="\\t"}}; {{print $1, $2 + {wildcards.cutbp}, $3 - {wildcards.cutbp} }}' {input.bed} > {output.bed}
        """

rule awk_cpgIslandExt_txt_to_bed3:
    """
    Created:
        2017-05-11 09:43:38
    Aim:

    Test:
        out/awk/cpgIslandExt_txt_to_bed3/gunzip/to-stdout/rsync/ucsc/goldenPath/mm10/database/cpgIslandExt.bed
    """
    input:
        "out/{filler}.txt"
    output:
        "out/awk/cpgIslandExt_txt_to_bed3/{filler}.bed"
    shell:
        """
        awk 'BEGIN{{FS=OFS="\\t"}} {{print $2,$3,$4}}' {input} > {output}
        """

rule awk_extract_bed6_from_rmsk:
    """
    Created:
        2018-08-16 20:12:15
    Aim:
        Exctracting bed6 file from rmsk
    Test:
        out/awk/extract_bed6_from_rmsk/gunzip/to-stdout/wget/http/hgdownload.cse.ucsc.edu/goldenpath/hg38/database/rmsk.bed

    """
    input:
        txt="out/{filler}.txt"
    output:
        bed="out/awk/extract_bed6_from_rmsk/{filler}.bed"
    threads:
        1
    shell:
        """
        awk 'BEGIN{{FS=OFS="\\t"}}{{
            {{print $6, $7, $8, $11 "__" $12 "__" $13, 0, $10}}
            }}' {input.txt} > {output.bed}
        """

rule awk_extract_bed3_from_rmsk_repeat_type:
    """
    Created:
        2017-01-09 11h56 - Adapted to new pattern concept and added wanted feature in output.
    Modified:
        2017-05-11 11:38:07
    Aim:
        This rule split the repeatMasker bed file into individual bed file for each repeat type.

    We select the 4th column in repeatMasker file which contain the name of the repeat type, then we keep only one of each type for the loop. Anchors tabulations are mandatory for grep because a lot of repeat type names are variations of the same pattern, e.g. searching L1ME3 will also retrieve regions called L1ME3A, L1ME3B and L1ME3C.

    Test:
        out/awk/extract_bed3_from_rmsk_repeat-type-GSAT_MM/gunzip/to-stdout/rsync/ucsc/goldenPath/mm10/database/rmsk.bed
    """
    input:
        txt="out/{filler}.txt"
    output:
        bed="out/awk/extract_bed3_from_rmsk_repeat-type-{feature}/{filler}.bed"
    threads:
        1
    wildcard_constraints:
        feature="[a-zA-Z0-9_-]+"
    shell:
        """
        awk 'BEGIN{{FS=OFS="\\t"}}{{
            if ($11 == "{wildcards.feature}"){{print $6,$7,$8}}
            }}' {input.txt} > {output.bed}
        """

rule awk_extract_bed6_from_cytoband:
    """
    Created:
        2017-07-27 16:06:15
    Aim:
        Cytoband txt file contains Giemsa stain results in 6th column. This rule produces a file in bed6 format converting gieStain value into score (gpos value * 6).
    Test:
        out/awk/extract_bed6_from_cytoband/gunzip/to-stdout/rsync/ucsc/goldenPath/mm10/database/cytoBand.bed
    Note:
        Possible values: gneg, gpos25, gpos33, gpos50, gpos66, gpos75, gpos100, acen, gvar, stalk
    """
    input:
        txt="out/{filler}.txt"
    output:
        bed="out/awk/extract_bed6_from_cytoband/{filler}.bed",
    threads:
        1
    shell:
        """
        awk 'BEGIN{{FS=OFS="\\t"}}{{
            gsub("gneg", "0", $5);
            gsub("gpos25", "150", $5);
            gsub("gpos33", "198", $5);
            gsub("gpos50", "300", $5);
            gsub("gpos66", "396", $5);
            gsub("gpos75", "450", $5);
            gsub("gpos100", "600", $5);
            gsub("acen", "750", $5);
            gsub("gvar", "875", $5);
            gsub("stalk", "1000", $5);
            print $1,$2,$3,$4,$5,"."
            }}' {input.txt} > {output.bed}
        """

rule awk_extract_from_bed4_feature_legacy:
    """
    Created:
        2017-06-19 15:48:45
    Modified:
        2018-11-13 16:18:38 - Changed to legacy because having the feature as basename is often more convenient.
    Aim:
        First used to prepare ChromHMM features to peak_anno
    Test:
        out/awk/extract_from_bed4_feature-E12/bedtools/makewindows_w-200_i-src/ChromHMM/LearnModel_test5_numstates-20_assembly-mm10/spermatozoa_20_segments.bed
    """
    input:
        txt="out/{filler}.bed"
    output:
        bed="out/awk/extract_from_bed4_feature-{feature}/{filler}.bed"
    threads:
        1
    wildcard_constraints:
        feature="[a-zA-Z0-9_-]+"
    shell:
        """
        awk 'BEGIN{{FS=OFS="\\t"}}{{
            if ($4 == "{wildcards.feature}"){{print $1,$2,$3}}
            }}' {input.txt} > {output.bed}
        """

rule awk_extract_from_bed4_feature:
    """
    Created:
        2017-06-19 15:48:45
    Aim:
        First used to prepare ChromHMM features to peak_anno
    Test:
        out/awk/extract_from_bed4_feature-E12/bedtools/makewindows_w-200_i-src/ChromHMM/LearnModel_test5_numstates-20_assembly-mm10/spermatozoa_20_segments.bed
    """
    input:
        txt="out/{filler}.bed"
    output:
        bed="out/awk/extract_feature_from_bed4/{filler}/{feature}.bed"
    threads:
        1
    wildcard_constraints:
        feature="[a-zA-Z0-9_-]+"
    shell:
        """
        awk 'BEGIN{{FS=OFS="\\t"}}{{
            if ($4 == "{wildcards.feature}"){{print $1,$2,$3}}
            }}' {input.txt} > {output.bed}
        """

rule awk_convert_bedpe_to_bed6_insert_size:
    """
    Created:
        2017-05-17 10:26:29
    Aim:
        This rule take bedpe file produced by bedtools bamtobed and convert it into a bed6 file where start and end correspond to the insert boundaries.
    Test:
        out/awk/convert_bedpe_to_bed6_insert_size/bedtools/bamtobed_bedpe/samtools/sort/samtools/merge/samtools/sort/samtools/view_bSh/bwa/mem_pe_GRCm38/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_en/inp/fastq/run163_run167/MNase_Spm_WT_band2.bed
    """
    input:
        bed="out/{filler}.bed"
    output:
        bed="out/awk/convert_bedpe_to_bed6_insert_size/{filler}.bed"
    threads:
        1
    shell:
        """
        awk 'BEGIN{{FS=OFS="\\t"}}{{
            if ($1 == $4 && $1!="."){{
                min=$2;
                if ($5 < $2){{
                    min=$5
                    }}
                max=$3
                if ($6 > $3){{
                    max=$6
                    }}
                print $1,min,max,$7,$8,"."
                }}
            }}' {input.bed} > {output.bed}
        """

#rule awk_select_subpopulation_from_sam_lmin_lmax:
#    input:
#        sam="out/{filler}.sam"
#    output:
#        sam="out/awk/select_subpopulations_from_sam_lmin-{lmin}_lmax-{lmax}/{filler}.sam"
#    wildcard_constraints:
#        lmin="[0-9]+",
#        lmax="[0-9]+"
#    shell:
#        """
#        awk 'BEGIN{{FS=OFS="\\t"}}{{
#            insert_length = sqrt($9^2)
#
#        """

rule awk_select_subpopulations_from_bed_lmin_lmax:
    """
    Created:
        2017-05-24 16:31:39
    Aim:
        Replace the selection of subpopulation made with my java code because it produces empty files in gtftk peak_anno.
    Note:
        out/awk/extract_main_chr/awk/convert_bedpe_to_bed6_insert_size/bedtools/bamtobed_bedpe/samtools/sort_n/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band3.bed
        out/awk/select_subpopulations_from_bed_lmin-100_lmax-130/awk/convert_bedpe_to_bed6_insert_size/bedtools/bamtobed_bedpe/samtools/sort_n/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band2.bed
        out/awk/select_subpopulations_from_bed_lmin-100_lmax-130/inp/unit_tests/awk_select_subpopulations_from_bed_lmin_lmax/1.bed

        out/awk/select_subpopulations_lmin-30_lmax-100/bowtie2/pe_mm10/ln/pe_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run163_run167_run184/s1-MNase_Spm_WT_band2_s2-MNase_Spm_WT_band4_s3-MNase_Spm_WT_band5.sam

        out/samtools/view_sam_to_bam/awk/select_subpopulations_lmin-30_lmax-100/bowtie2/pe_hg19/sickle/pe_-t_sanger_-q_20/gunzip/already_merged_lane_nextseq500_pe/ln/updir/mw/inp/fastq/run246/S002740_803_K4me3-113115.bam
    """
    input:
        bed="out/{filler}.{ext}"
    output:
        bed="out/awk/select_subpopulations_lmin-{lmin}_lmax-{lmax}/{filler}.{ext}"
    params:
        insert_length = lambda wildcards: "$3 - $2" if wildcards.ext == 'bed' else "sqrt($9^2)"
    wildcard_constraints:
        ext="bed|sam",
        lmin="[0-9]+",
        lmax="[0-9]+"
    shell:
        """
        awk 'BEGIN{{FS=OFS="\\t"}}{{
            insert_length = {params.insert_length}
            if ($0 ~ /^@/){{
                print $0
            }}
            if (insert_length < {wildcards.lmax} && insert_length > {wildcards.lmin}){{
                print $0
                }}
            }}' {input.bed} > {output.bed}
        """

rule awk_select_subpopulations_from_bed_lmin:
    """
    Created:
        2018-11-23 14:58:43
    Aim:
    """
    input:
        bed="out/{filler}.bed"
    output:
        bed="out/awk/select_subpopulations_from_bed_lmin-{lmin}/{filler}.bed"
    wildcard_constraints:
        lmin="[0-9]+",
    shell:
        """
        awk 'BEGIN{{FS=OFS="\\t"}}{{
            insert_length = $3 - $2
            if (insert_length > {wildcards.lmin}){{
                print $0
                }}
            }}' {input.bed} > {output.bed}
        """

rule awk_select_subpopulations_from_bed_lmax:
    """
    Created:
        2018-11-23 14:58:43
    Aim:
    """
    input:
        bed="out/{filler}.bed"
    output:
        bed="out/awk/select_subpopulations_from_bed_lmax-{lmax}/{filler}.bed"
    wildcard_constraints:
        lmax="[0-9]+",
    shell:
        """
        awk 'BEGIN{{FS=OFS="\\t"}}{{
            insert_length = $3 - $2
            if (insert_length < {wildcards.lmax}){{
                print $0
                }}
            }}' {input.bed} > {output.bed}
        """

rule awk_extract_xls_coordinates_to_bed3:
    """
    Created:
        2017-01-24 17h09
    Test:
        out/awk/extract_xls_coordinates_to_bed3/macs2/callpeak_broad_BAMPE_mm_SPMRTrue/bam/mm9/merge_run140_run141/H4K5ac-Nut-WT_over_Input-Nut-WT_peaks.bed
        out/awk/extract_xls_coordinates_to_bed3/macs2/callpeak_broad_BAM_mm_SPMRTrue/samtools/sam_to_bam/bowtie2/se/mm9/GSE56526/SRR1596612_over_SRR1596620_peaks.bed
        out/awk/extract_xls_coordinates_to_bed3/macs2/callpeak_broad_qvalue0.01_BAM_mm/samtools/sam_to_bam/bowtie2/se/mm9/GSE56526/SRR1596612_over_SRR1596620_peaks.bed
    """
    input:
        xls="out/{filler}.xls"
    output:
        bed="out/awk/extract_xls_coordinates_to_bed3/{filler}.bed"
    shell:
        """
        awk 'BEGIN{{FS=OFS="\\t"}} $1 ~ /^chr([0-9]+|X|Y)$/ {{print $1,$2,$3}}' {input.xls} > {output.bed}
        """

rule awk_convert_sicer_to_bed6:
    """
    Created:
        2018-03-12 18:10:40
    Aim:
        Because SICER output is not true bed, I need this conversion to then convert to GFF.
    Test:
        out/awk/convert_sicer_to_bed6/awk/extract_main_chr/crossmap/bed_mm9_to_mm10/gunzip/tar/xvf/data_brdt/GSM984200_4_EM2_R_BRDT-Paired.F3.no_4-W200-G600-islands-summary-FDR001.bed out/awk/extract_main_chr/crossmap/bed_mm9_to_mm10/awk/extract_xls_coordinates_to_bed3/macs2/callpeak_broad_BAM_mm_SPMRTrue/samtools/sam_to_bam/bowtie2/se/mm9/GSE56526/SRR1596612_over_SRR1596620_peaks.bed
    """
    input:
        bed="out/{filler}.bed"
    output:
        bed="out/awk/convert_sicer_to_bed6/{filler}.bed"
    conda:
        "../envs/gawk.yaml"
    shell:
        """
        awk 'BEGIN{{FS=OFS="\\t"}} {{print $1,$2,$3,".",".","+"}}' {input.bed} > {output.bed}
        """

rule awk_convert_bed3_to_saf:
    """
    Created:
        2019-03-03 20:42:25
    Aim:
        Convert bed3 to saf using the merge of chr, start and end as new feature name.
    Test:
        out/awk/convert_bed3_to_saf/macs2/callpeak_--broad/samtools/index/samtools/sort/samtools/view_sam_to_bam_-q_30/bowtie2/se_mm10/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR3126243_over_SRR3126242_peaks.saf
    """
    input:
        bed="out/{filler}.bed"
    output:
        saf="out/awk/convert_bed3_to_saf/{filler}.saf"
    conda:
        "../envs/gawk.yaml"
    shell:
        """
        awk 'OFS="\t" {{print $1"-"$2+1"-"$3, $1, $2+1, $3, "+"}}' {input} > {output.saf}
        """

rule awk_convert_bed6_to_saf:
    """
    Created:
        2019-03-03 20:42:25
    Aim:
        Convert bed6 to saf
    Test:
        out/awk/convert_bed6_to_saf/bedtools/slop_-b_1000_chrominfo-mm10/gtftk/get_5p_3p_coords/gtftk/select_most_5p_tx/awk/extract_main_chr/sed/add_chr/gunzip/to-stdout/wget/ftp/ftp.ensembl.org/pub/release-95/gtf/mus_musculus/Mus_musculus.GRCm38.95.bed
    """
    input:
        bed="out/{filler}.bed"
    output:
        saf="out/awk/convert_bed6_to_saf/{filler}.saf"
    conda:
        "../envs/gawk.yaml"
    shell:
        """
        awk 'OFS="\t" {{print $4, $1, $2+1, $3, $6}}' {input} > {output.saf}
        """




rule awk_convert_macs2_peaks_to_great:
    """
    Created:
        2018-03-13 10:27:22
    Aim:
        Because GREAT is picky on bed format
    Note:
         Line 1: The seventh field (thickStart) must be a non-negative integer (not 3.12785). The eighth field (thickEnd) must be a non-negative integer (not 4.07768). The ninth field (itemRgb) must be either a non-negative integer or a comma-delimited list of 3 non-negative integers (not 1.26763). Input data: 'chr1 29335341 29336038 Th91_SP4_over_input_peak_203 12 . 3.12785 4.07768 1.26763'

         chr1    29335341        29336038        Th91_SP4_over_input_peak_203    12      .       3.12785 4.07768 1.26763

    """
    input:
        bed="out/{filler}.bed"
    output:
        bed="out/awk/convert_macs2_peaks_to_great/{filler}.bed"
    conda:
        "../envs/gawk.yaml"
    shell:
        """
        awk 'BEGIN{{FS=OFS="\\t"}} {{print $1,$2,$3,$4,$5,$6}}' {input.bed} > {output.bed}
        """

rule awk_convert_to_bed6:
    """
    Created:
        2018-03-13 22:58:53
    Aim:
        Because it seems crossmap does not work on bed12, so I try to just use bed6.

    """
    input:
        bed="out/{filler}.bed"
    output:
        bed="out/awk/convert_to_bed6/{filler}.bed"
    conda:
        "../envs/gawk.yaml"
    shell:
        """
        awk 'BEGIN{{FS=OFS="\\t"}} {{print $1,$2,$3,$4,"0",$6}}' {input.bed} > {output.bed}
        """

rule awk_extract_macs2_xls_coordinates_to_bed6:
    """
    Created:
        2017-02-21 17:11:49
    Aim:
        Because bed3 is not a true bed and can not be taken by some tools...
        $8 corresponds to score column in macs2 peaks but should not work with other tools.
    Test:
        out/awk/extract_xls_coordinates_to_bed6/macs2/callpeak_broad_BAMPE_mm_SPMRTrue/bam/mm9/merge_run140_run141/H4K5ac-Nut-WT_over_Input-Nut-WT_peaks.bed
        out/awk/extract_xls_coordinates_to_bed6/macs2/callpeak_broad_BAM_mm_SPMRTrue/samtools/sam_to_bam/bowtie2/se/mm9/GSE56526/SRR1596612_over_SRR1596620_peaks.bed
        out/awk/extract_xls_coordinates_to_bed6/macs2/callpeak_broad_qvalue0.01_BAM_mm/samtools/sam_to_bam/bowtie2/se/mm9/GSE56526/SRR1596612_over_SRR1596620_peaks.bed
        out/awk/extract_xls_coordinates_to_bed6/macs2/callpeak_broad_format-BAM_gsize-hs_keepdup-all_no_control/ln/alias/experiments/hg38_ATAC_thymus/CD34_peaks.bed
    """
    input:
        xls="out/{filler}.xls"
    output:
        bed="out/awk/extract_xls_coordinates_to_bed6/{filler}.bed"
    shell:
        """
        awk 'BEGIN{{FS=OFS="\\t"}} $1 ~ /^chr([0-9]+|X|Y)$/ {{print $1,$2,$3,NR,$8,"."}}' {input.xls} > {output.bed}
        """

rule awk_extract_genes_from_gtftk_convert:
    """
    Created:
        2017-09-27 17:23:12
    Modified:
        2018-11-01 09:53:45 - The pattern has changed in recente pygtftk version, so '|\?' is used instead of '|\.'.
    Aim:
        When converting gtf to bed, "gtftk convert" produces one line per gene and one line for each transcript. If this bed file has to be given to deepTools, it is better to keep only the line with the gene coordinate. Such lines have their bed name ending with "|." which is convenient.

    Test:
        out/awk/extract_genes_from_gtftk_convert/gtftk/convert/gtftk/select_by_key_key-gene_id_file-with-values-upreg-ko-nut-r-threshold-1/gunzip/to-stdout/wget/ftp_ensembl/pub/release-87/gtf/mus_musculus/Mus_musculus.GRCm38.87.bed
    """
    input:
        bed="out/{filler}.bed"
    output:
        bed="out/awk/extract_genes_from_gtftk_convert/{filler}.bed"
    shell:
        """
        awk 'BEGIN{{FS=OFS="\\t"}} $4 ~ /^*|\?$/ {{print $0}}' {input.bed} | sed 's/|\.//g' > {output.bed}
        """

rule awk_extract_gfold_signif_threshold:
    """
    Created:
        2017-08-02 22:16:40
    Modified:
        2017-09-27 16:13:10 - Renamed from awk_extract_gfold_signif_peaks to awk_extract_gfold_signif_threshold.
    Aim:
        Extract significant positions to bed in order to use them as reference for heatmap.
    Test:
        input:
            out/r/sort_gfold_diff_by_gfold_then_log2ratio/gfold/diff/gfold/count_bed-mm10-R-H4K5ac-bu-peaks/samtools/sort/samtools/view_bSh/bowtie2/se_mm10/sickle/se_-t_sanger_-q_20/sra-tools/fastqdump_se/wget/ftp_trace_ncbi_sra/SRR3126243_VS_SRR3126244_diff_ascend.tsv
        output:
            out/awk/extract_gfold_signif_threshold-0/gfold/diff/gfold/count_bed-mm10-R-H4K5ac-bu-peaks/samtools/sort/samtools/view_bSh/bowtie2/se_mm10/sickle/se_-t_sanger_-q_20/sra-tools/fastqdump_se/wget/ftp_trace_ncbi_sra/SRR3126243_VS_SRR3126244_pos.bed

            out/awk/extract_gfold_signif_threshold-1/gfold/diff/gfold/count_gtf-GRCm38/star/pe_GRCm38/sickle/pe_-t_sanger_-q_20/gunzip/nut_rnaseq/Nut-R-WT_VS_Nut-R-KO_pos.bed
    """
    input:
        diff="out/r/sort_gfold_diff_by_gfold_then_log2ratio/{filler}_diff_ascend.tsv"
    output:
        txt_gfold_pos="out/awk/extract_gfold_signif_threshold-{threshold}/{filler}_pos.txt",
        txt_gfold_neg="out/awk/extract_gfold_signif_threshold-{threshold}/{filler}_neg.txt"
    wildcard_constraints:
        threshold='0|1|2|3|4|5'
    shell:
        """
        awk 'BEGIN {{OFS="\\t"}}; {{if ($3>{wildcards.threshold}) print $1}}' {input.diff} | sed 's/"//g' | sed 's/-/\\t/g' > {output.txt_gfold_pos}
        awk 'BEGIN {{OFS="\\t"}}; {{if ($3<-{wildcards.threshold}) print $1}}' {input.diff} | sed 's/"//g'| sed 's/-/\\t/g' > {output.txt_gfold_neg}
        """

rule awk_extract_nut_microarray_signif_threshold:
    """
    Created:
        2017-10-02 14:56:53
    Aim:
        Extracts genes ids from sophie rousseaux shared results from microarrays.
    Test:
        input:
            out/sed/mac_to_unix/inp/sophie_rousseaux_2017_09_29/anova_nut_sp_line_res_3_RSkovRSwt_gsea1.txt
        output:
            out/awk/extract_nut_microarray_signif_threshold-1.5/sed/mac_to_unix/inp/sophie_rousseaux_2017_09_29/anova_nut_sp_line_res_3_RSkovRSwt_gsea1_pos.txt
    """
    input:
        txt="out/{filler}.txt"
    output:
        txt_signif_pos="out/awk/extract_nut_microarray_signif_threshold-{threshold}/{filler}_pos.txt",
        txt_signif_neg="out/awk/extract_nut_microarray_signif_threshold-{threshold}/{filler}_neg.txt"
    wildcard_constraints:
        threshold='1.5|1.2'
    shell:
        """
        awk 'BEGIN {{FS=OFS="\\t"}}; {{if ($2>{wildcards.threshold}) print $1}}' {input.txt}  > {output.txt_signif_pos}
        awk 'BEGIN {{FS=OFS="\\t"}}; {{if ($2<-{wildcards.threshold}) print $1}}' {input.txt} > {output.txt_signif_neg}
        """

rule awk_extract_chr9:
    """
    Created:
        2017-11-02 14:04:04
    Aim:
        Extract lines starting with chr9 for test.
    Test:
    """
    input:
        bed="out/{filler}.bed"
    output:
        bed="out/awk/extract_chr9/{filler}.bed"
    shell:
        """
        awk 'BEGIN {{FS=OFS="\\t"}}; {{if ($1 == "chr9") print $0}}' {input.bed}  > {output.bed}
        """

rule awk_extract_genes_from_gff:
    """
    Created:
        2018-03-08 11:41:56
    Aim:
        Extract lines starting with chr9 for test.
    Test:
        out/awk/extract_genes_from_gff/gunzip/to-stdout/wget/ftp/ftp.ensembl.org/pub/release-91/gff3/mus_musculus/Mus_musculus.GRCm38.91.chromosome.19.gff3
    """
    input:
        gff="out/{filler}"
    output:
        gff="out/awk/extract_genes_from_gff/{filler}"
    shell:
        """
        awk 'BEGIN {{FS=OFS="\\t"}}; {{if ($3 == "gene") print $0}}' {input.gff}  > {output.gff}
        """

"""
chr1    gffFromBW       gffFromBW       901     1001    0       +       .       .
chr19   ensembl_havana  gene    3260924 3283017 .       -       .       ID=gene:ENSMUSG00000024831;Name=Ighmbp2;biotype=protein_coding;description=immunoglobulin mu binding protein 2 [Source:MGI Symbol%3BAcc:MGI:99954];gene_id=ENSMUSG00000024831;havana_gene=OTTMUSG00000028273;havana_version=1;logic_name=ensembl_havana_gene;version=13

"""

rule awk_make_gff_looks_like_example_from_chipseqspike:
    """
    Created:
        2018-03-08 14:10:04
    Test:
        out/awk/make_gff_looks_like_example_from_chipseqspike/sed/add_chr/awk/extract_genes_from_gff/gunzip/to-stdout/wget/ftp/ftp.ensembl.org/pub/release-91/gff3/mus_musculus/Mus_musculus.GRCm38.91.chromosome.19.gff
    """
    input:
        gff="out/{filler}.gff3"
    output:
        gff="out/awk/make_gff_looks_like_example_from_chipseqspike/{filler}.gff"
    shell:
        """
        awk 'BEGIN {{FS=OFS="\\t"}}; {{print $1,$3,$3,$4,$5,0,$7,".","."}}' {input.gff} > {output.gff}
        """

rule awk_extract_rseg_score_above_threshold:
    """
    Created:
        2017-11-03 11:26:24
    Aim:
        Extract differential peaks with score above threshold because sometimes the peaks are questionable.
    Note:
        After testing, it does not seem that score is the best discriminant for "good" differential peaks.
    Test:
        out/awk/extract_rseg_score_above_threshold-10/rseg_diff_chrom-mm9_deadzones-mm9_mode-3/sort/coordinates_rseg_bed/awk/convert_bedpe_to_bed6_insert_size/bedtools/bamtobed_bedpe/samtools/sort_n/samtools/merge/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_mm9/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run140_run141/H4K5ac-Nut-WT_VS_H4K5ac-Nut-KO.bed
    """
    input:
        bed="out/{filler}.bed"
    output:
        bed="out/awk/extract_rseg_score_above_threshold-{threshold}/{filler}.bed"
    wildcard_constraints:
        threshold="[0-9]+"
    shell:
        """
        awk 'BEGIN {{FS=OFS="\\t"}}; {{if ($6>{wildcards.threshold}) print $0}}' {input.bed}  > {output.bed}
        """

rule awk_extract_rseg_avg_count_diff_above_threshold:
    """
    Created:
        2017-11-03 14:00:26
    Aim:
        Extract differential peaks with average count difference above threshold because sometimes the peaks are questionable. Also the score alone does not extract the most humanly trustable peaks.
    Test:
        out/awk/extract_rseg_avg_count_diff_above_threshold-20/rseg_diff_chrom-mm9_deadzones-mm9_mode-3/sort/coordinates_rseg_bed/awk/convert_bedpe_to_bed6_insert_size/bedtools/bamtobed_bedpe/samtools/sort_n/samtools/merge/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_mm9/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run140_run141/H4K5ac-Nut-WT_VS_H4K5ac-Nut-KO_pos.bed
    """
    input:
        txt="out/{filler}.bed"
    output:
        bed_pos="out/awk/extract_rseg_avg_count_diff_above_threshold-{threshold}/{filler}_pos.bed",
        bed_neg="out/awk/extract_rseg_avg_count_diff_above_threshold-{threshold}/{filler}_neg.bed"
    wildcard_constraints:
        threshold="[0-9]+"
    shell:
        """
        awk 'BEGIN {{FS=OFS="\\t"}}; {{if ($5> {wildcards.threshold}) print $0}}' {input.txt} > {output.bed_pos}
        awk 'BEGIN {{FS=OFS="\\t"}}; {{if ($5<-{wildcards.threshold}) print $0}}' {input.txt} > {output.bed_neg}
        """

rule awk_extract_rseg_avg_count_diff_and_domain_score_above_threshold:
    """
    Created:
        2017-11-03 14:00:26
    Aim:
        Extract differential peaks with average count difference above threshold because sometimes the peaks are questionable. Also the score alone does not extract the most humanly trustable peaks.
    Test:
        out/awk/extract_rseg_avg_count_diff_and_domain_score_above_threshold-10/rseg_diff_chrom-mm9_deadzones-mm9_mode-3/sort/coordinates_rseg_bed/awk/convert_bedpe_to_bed6_insert_size/bedtools/bamtobed_bedpe/samtools/sort_n/samtools/merge/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_mm9/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run140_run141/H4K5ac-Nut-WT_VS_H4K5ac-Nut-KO_pos.bed
    """
    input:
        txt="out/{filler}.bed"
    output:
        bed_pos="out/awk/extract_rseg_avg_count_diff_and_domain_score_above_threshold-{threshold}/{filler}_pos.bed",
        bed_neg="out/awk/extract_rseg_avg_count_diff_and_domain_score_above_threshold-{threshold}/{filler}_neg.bed"
    wildcard_constraints:
        threshold="[0-9]+"
    shell:
        """
        awk 'BEGIN {{FS=OFS="\\t"}}; {{if ($5> {wildcards.threshold} && $6>{wildcards.threshold}) print $0}}' {input.txt} > {output.bed_pos}
        awk 'BEGIN {{FS=OFS="\\t"}}; {{if ($5<-{wildcards.threshold} && $6>{wildcards.threshold}) print $0}}' {input.txt} > {output.bed_neg}
        """

rule awk_rseg_domains_to_bed:
    """
    Created:
        2017-11-03 14:39:16
    Aim:
        Because RSEG domains have two score columns (avg_count_diff_and_domain_score) and this is not cool for tools like GREAT. So here I merge these two scores with multiplication and put the result as integer.
    Test:
        out/awk/rseg_domains_to_bed/awk/extract_rseg_avg_count_diff_and_domain_score_above_threshold-10/rseg_diff_chrom-mm9_deadzones-mm9_mode-3/sort/coordinates_rseg_bed/awk/convert_bedpe_to_bed6_insert_size/bedtools/bamtobed_bedpe/samtools/sort_n/samtools/merge/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_mm9/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run140_run141/H4K5ac-Nut-WT_VS_H4K5ac-Nut-KO_pos.bed
        out/awk/rseg_domains_to_bed/awk/extract_rseg_avg_count_diff_and_domain_score_above_threshold-1/rseg_diff_chrom-mm9_deadzones-mm9_mode-3/sort/coordinates_rseg_bed/awk/convert_bedpe_to_bed6_insert_size/bedtools/bamtobed_bedpe/samtools/sort_n/samtools/merge/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_mm9/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run140_run141/H4K5ac-Nut-WT_VS_H4K5ac-Nut-KO_pos.bed
    """
    input:
        bed="out/{filler}.bed"
    output:
        bed="out/awk/rseg_domains_to_bed/{filler}.bed"
    wildcard_constraints:
        threshold="[0-9]+"
    shell:
        """
        awk 'BEGIN {{FS=OFS="\\t"}}; {{print $1,$2,$3,$4,int($5*$6),$7}}' {input.bed}  > {output.bed}
        """

rule awk_extract_distal_proximal:
    """
    Created:
        2018-03-13 18:32:17
    Note:
        ==> out/paste/great_input_and_output_hg19/crossmap/bed_hg38_to_hg19/bedtools/intersect_thymus_stage_common_peaks_hg38_v3/not_common_subset_3.bed <==
        chr1    713374  714740  92004
        chr1    762103  763261  98436
        chr1    773467  774020  87374
        chr1    776556  777424  84128
        chr1    778226  780332  81839
        chr1    839685  840550  21000
        chr1    845517  849555  13582
        chr1    877416  877757  16469
        chr1    893948  896689  648
        chr1    901745  902363  177
    Test:
        out/awk/extract_distal_proximal/out/paste/great_input_and_output_hg19/crossmap/bed_hg38_to_hg19/bedtools/intersect_thymus_stage_common_peaks_hg38_v3/not_common_subset_3.bed
    """
    input:
        bed="out/{filler}.bed"
    output:
        distal="out/awk/extract_distal_proximal/{filler}_distal.bed",
        proximal="out/awk/extract_distal_proximal/{filler}_proximal.bed",
    shell:
        """
        awk 'BEGIN {{FS=OFS="\\t"}}; {{if ($NF > 2000) print $0}}' {input.bed} | cut -f1-6 > {output.distal}
        # Additionnal check if last field is not empty, because empty field count as number below 2000 in awk...
        awk 'BEGIN {{FS=OFS="\\t"}}; {{if ($NF < 2000 && $NF != "") print $0}}' {input.bed} | cut -f1-6 > {output.proximal}
        """

#rule awk_bedtools_closest_to_bed6:
#    """
#    Created:
#        2018-09-05 23:54:28
#    Aim:
#        Convert this to bed6 to import in R with rtracklayer
#        $ head out/bedtools/closest_d_b-hg38-ensembl-r93-tss-protein-coding-TR-IG/sort/coordinates_bed/bedtools/merge/bedtools/multiinter_bed-hg38-bs-hypometh-thymus/multiinter.bed
#        chr1    29337   29393   chr1    65418   65419   ENSG00000186092|ENST00000641515 .       +       36026
#        chr1    91267   91581   chr1    69054   69055   ENSG00000186092|ENST00000335137 .       +       22213
#        chr1    191747  191777  chr1    69054   69055   ENSG00000186092|ENST00000335137 .       +       122693
#
#    """

rule awk_extract_plotHeatmap_cluster:
    """
    Created:
        2018-03-13 22:16:10
    Aim:
        After a k-means done by deepTools plotHeatmap, I want to extract bed containing only one cluster to be able to submit this cluster to GREAT analysis.
    Note:
        cluster is in column 13.
    Test:
        input:
            out/deepTools/plotHeatmap_kmeans-10_sortRegions-descend_sortUsing-region_length_averageTypeSummaryPlot-mean_colorList-blueCyanYellowOrangeRed_heatmapHeight-28_heatmapWidth-4_whatToShow-phc_xAxisLabel-tss_refPointLabel-0/deepTools/computeMatrix_reference-point_referencePoint-center_beforeRegionStartLength-2000_afterRegionStartLength-2000_bed-hg38-macs2-peaks-H3K27ac-thymus-stage-exclusive-not-common-subset-3-distal_bw-hg38-H3K27ac-thymus-subset-3.bed
        output:
            out/awk/extract_plotHeatmap_cluster/deepTools/plotHeatmap_kmeans-10_sortRegions-descend_sortUsing-region_length_averageTypeSummaryPlot-mean_colorList-blueCyanYellowOrangeRed_heatmapHeight-28_heatmapWidth-4_whatToShow-phc_xAxisLabel-tss_refPointLabel-0/deepTools/computeMatrix_reference-point_referencePoint-center_beforeRegionStartLength-2000_afterRegionStartLength-2000_bed-hg38-macs2-peaks-H3K27ac-thymus-stage-exclusive-not-common-subset-3-distal_bw-hg38-H3K27ac-thymus-subset-3_c-9.bed
    """
    input:
        bed="out/{filler}.bed"
    output:
        bed="out/awk/extract_plotHeatmap_cluster/{filler}_c-{cluster}.bed"
    wildcard_constraints:
        cluster="[0-9]+"
    shell:
        """
        awk 'BEGIN {{FS=OFS="\\t"}}; {{if ($13 == "cluster_{wildcards.cluster}") print $0}}' {input.bed} > {output.bed}
        """

#rule awk_extract_cluster_to_bed6:
#    """
#    Created:
#        2018-02-22 14:25:12
#    Aim:
#    Test:
#    """
#    input:
#        bed="out/{filler}.bed"
#    output:
#        bed="out/awk/extract_cluster-{cluster}_to_bed6/{filler}.bed"
#    wildcard_constraints:
#        cluster="[0-9]+"
#    shell:
#        """
#        awk 'BEGIN {{OFS="\\t"}}; ($5 > {wildcards.minsmt}) && ($6 < {wildcards.maxfuz}) {{print $1, $2, $3;}}' {input.xls} > {output.bed}
#        """
#
#
rule awk_filter_bedgraph_gt_threshold:
    """
    Created: 2016-12-20 23h08 - Trying to see coverage stats we could get if instead of threshold 1 we used higher threshold.

    Test:
        "out/bedtools/genomecov/bedgraph/awk/filter_bedgraph/gt3/bedtools/genomecov/bed_to_bedgraph/awk/filter_bedpe_by_size/min40_max80/samtools/view/bampe_to_bed3_nodup/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126/PSK-SC-KO.txt"
    """
    input:
        bedgraph="out/{filler}.bedgraph"
    output:
        bedgraph="out/awk/filter_bedgraph/gt{threshold}/{filler}.bedgraph"
    wildcard_constraints:
        threshold="[0-9]+"
    shell:
        """
        awk 'BEGIN {{OFS="\\t"}}; {{if ($4 > {wildcards.threshold}) print $0}}' {input.bedgraph} > {output.bedgraph}
        """

rule awk_prepare_include_for_shuffle:
    """
    Created: 2016-12-07 17h02 - Here I use a mappability track retrieved from UCSC in bw format then converted in bedgraph to define the regions of the genome that are perfectly mappable (the 4th column contains '1'). This is aimed to produce more accurate shuffling for negative control.
    Maybe it would be more accurate to select regions with mappability greater than X instead of only 1, but I have no element to favour one approach right now.
    """
    input:
        "out/ucsc/bigWigToBedGraph/rsync/ucsc/gbdb/mm9/bbi/crgMapabilityAlign75mer.bedGraph"
    output:
        "out/awk/prepare_include_for_shuffle/ucsc/bigWigToBedGraph/rsync/ucsc/gbdb/mm9/bbi/crgMapabilityAlign75mer.bedGraph"
    shell:
        """
        awk 'BEGIN {{FS="\\t"}}; {{if ($4==1) print $0}}' {input} > {output}
        """

rule awk_get_non_unmapable:
    """
    Created:
        2016-12-16 14h18 - Because using include does not seem to work very well with bedtools shuffle, I try the opposite approach with exclude.
    """
    input:
        mapability="out/ucsc/bigWigToBedGraph/rsync/ucsc/gbdb/mm9/bbi/crgMapabilityAlign75mer.bedGraph"
    output:
        non_unmapable="out/awk/get_non_unmapable/ucsc/bigWigToBedGraph/rsync/ucsc/gbdb/mm9/bbi/crgMapabilityAlign75mer.bedGraph"
    shell:
        """
        awk 'BEGIN {{FS="\\t"}}; {{if ($4!=0) print $0}}' {input} > {output}
        """

rule awk_filter_main_chr:
    input:
        bedlike="out/{filler}"
    output:
        bedlike="out/awk/filter_main_chr/{filler}"
    shell:
        """
        awk '!/_random/' {input.bedlike} > {output.bedlike}
        """

#rule awk_chromInfo_to_bed:
#    input:
#        chromInfo="out/mysql/ucsc/mm9/chromInfo_main.tsv"
#    output:
#        bed="out/awk/chromInfo_to_bed/mysql/ucsc/mm9/chromInfo_main.bed"
#    shell:"""
#    awk 'BEGIN {{OFS="\\t"}}; {{print $1, 0, $2}}' {input.chromInfo} > {output.bed}
#    """

rule awk_filter_bedpe_by_size:
    """
    Created: 2016-12-20 15h52

    Test:
        "out/awk/filter_bedpe_by_size/min40_max80/samtools/view/bampe_to_bed3_nodup/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126/PSK-SC-KO.bed"
        "out/awk/filter_bedpe_by_size/min40_max80/samtools/view/bampe_to_bed3_nodup/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126/PSK-SC-KO.bed"


    """
    input:
        bed="out/{filler}.bed"
    output:
        bed="out/awk/filter_bedpe_by_size/min{min}_max{max}/{filler}.bed"
    wildcard_constraints:
        min="[0-9]+",
        max="[0-9]+"
    shell:
        """
        awk 'BEGIN {{OFS="\\t"}}; {{if (($3-$2 < {wildcards.max}) && ($3-$2 > {wildcards.min})) {{print $0}} }}' {input.bed} > {output.bed}
        """
