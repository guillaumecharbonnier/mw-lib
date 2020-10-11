"""
Note:
    Depending on your needs, you may rather use on of these three rules:
    bedtools_intersect_b_lambda_extra
    bedtools_intersect_a_lambda_extra
    bedtools_intersect_samedir_extra
Doc:
    https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html
"""

rule bedtools_intersect_b_lambda_extra:
    """
    Created:
        2017-11-14 17:21:40
    Test:
        out/bedtools/intersect_-wa_-b_bed-hg19-active-enhancers-thymopoiesis-tall-samples/sort/_-k1,1_-k2,2n/cat/cat-hg19-active-enhancers-thymopoiesis-tall-samples.bed
        out/bedtools/intersect_-b_bed-hg19-polycomb-in-at-least-4-thymocytes/r/active_tss_dynamics_hsc-tcell-tall-samples/at_least_3_TALL_no_thymocyte_tss.bed
    """
    input:
        features_a = "out/{filler}.{ext}",
        features_b = lambda wildcards: eval(config['ids'][wildcards.feature_list_id])
    output:
        features_a_overlapping_features_b="out/{tool}{extra}_-b_{feature_list_id}/{filler}.{ext}"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="bedtools/intersect",
        feature_list_id = "[a-zA-Z0-9-]+",
        ext="bam|bed|bedgraph|gff|vcf"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.features_a} -b {input.features_b} {params.extra} > {output.features_a_overlapping_features_b}"

rule bedtools_intersect_a_lambda_extra:
    """
    Created:
        2020-06-08 20:12:50
    Test:
    """
    input:
        features_b = "out/{filler}.{ext}",
        features_a = lambda wildcards: eval(config['ids'][wildcards.bed_list_id])
    output:
        features_a_overlapping_features_b="out/{tool}{extra}_-a_{bed_list_id}/{filler}.{ext}"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="bedtools/intersect",
        ext="bam|bed|bedgraph|gff|vcf"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.features_a} -b {input.features_b} {params.extra} > {output.features_a_overlapping_features_b}"

rule bedtools_intersect_samedir_extra:
    """
    Created:
        2020-06-06 23:47:32
    Aim:
        More convenient than rule above when both a and b are from the same directory.
    Test:
    """
    input:
        features_a = "out/{filler}/{a}.{ext}",
        features_b = "out/{filler}/{b}.{ext}"
    output:
        features_a_overlapping_features_b="out/{tool}{extra}/{filler}/{a}_VS_{b}.{ext}"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="bedtools/intersect",
        ext="bam|bed|bedgraph|gff|vcf"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.features_a} -b {input.features_b} {params.extra} > {output.features_a_overlapping_features_b}"

### Legacy past this point. Or very specific rules.

rule bedtools_intersect_wa:
    """
    Created:
        2017-11-14 17:21:40
    Test:
        features_a:
            out/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_GRCm38/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run124/MNS-R-WT.bam
        features_b:
            out/sh/danpos_xls_to_bed3/danpos/dtriple_v4/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run113_run119_run124/MNS-R-WT/danpos.smooth.positions.bed
        features_a_overlapping_features_b:
            out/bedtools/intersect_wa/b-sh/danpos_xls_to_bed3/danpos/dtriple_v4/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_GRCm38/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run113_run119_run124/MNS-R-WT/danpos.smooth.positions.bed/a-samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_GRCm38/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run124/MNS-R-WT.bam
        next step:
            out/r/dinucleotide_plot_u-300_d-300_f-GRCm38/perl/peReadsFragmentMid_f-147-147/bedtools/intersect_wa/b-sh/danpos_xls_to_bed3/danpos/dtriple_v4/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_GRCm38/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run113_run119_run124/MNS-R-WT/danpos.smooth.positions.bed/a-samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_GRCm38/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run124/MNS-R-WT.pdf

        out/bedtools/intersect_wa/b-out/crossmap/bed_hg19_to_hg38/awk/extract_distal_proximal/paste/great_input_and_output_hg19/crossmap/bed_hg38_to_hg19/awk/fill_bed3_to_bed6/bedtools/multiinter_thymus_peaks_hg38_merged/1_distal.bed/a-out/bedtools/multiinter_bed-hg38-bs-hypometh-thymus/common_all.bed
    """
    input:
        features_a="out/{filler1}.{ext1}",
        features_b="out/{filler2}.{ext2}"
    output:
        features_a_overlapping_features_b="out/bedtools/intersect_wa/b-{filler2}.{ext2}/a-{filler1}.{ext1}"
    wildcard_constraints:
        ext1="bam|bed|bedgraph|gff|vcf",
        ext2="bam|bed|bedgraph|gff|vcf"
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        bedtools\
            intersect\
            -a {input.features_a}\
            -b {input.features_b}\
            -wa > {output.features_a_overlapping_features_b}
        """



#def input_bed_for_rule_bedtools_intersect_filter_bam(wildcards):
#	feature_path= {'mm9'+'tss_pm2kb': 'mm9TSSRegions'}
#	out_path=wildcards['index']+wildcards['feature']
#        print(out_path)
#	return "annotation/input/"+ feature_path[out_path] +".bed"
#

#
#def input_bed_for_rule_bedtools_intersect_filter_bam(index='mm9',feature='tss_pm2kb'):
#        feature_path= {'mm9'+'tss_pm2kb': 'mm9TSSRegions'}
#        out_path=index+feature
#        return "annotation/input/"+ feature_path[out_path] +".bed"

rule bedtools_intersect_bed_exclude_feature:
    """
    Created: 2016-07-21 - Made to refine the features used to spot differences in fragment size distribution between small structures coming from nucleosomes or from TSS regions (need manual preparation to look for NFR only)

    Usage:
    "annotation/processed/feature/mm9/danpos/run119/MNS-R-WT/positions_exclude_tss_pm2kb.bed",
    "annotation/processed/feature/mm9/tss_pm2kb_exclude_danpos/run119/MNS-R-WT/positions.bed",

    """
    input:
        bed_coordinates="annotation/processed/feature/{index}/{feature}.bed",
        bed_feature_to_exclude="annotation/processed/feature/{index}/{feature_to_exclude}.bed"
    output:
        bed="annotation/processed/feature/{index}/{feature, danpos/run119/MNS-R-WT/positions|tss_pm2kb}_exclude_{feature_to_exclude, danpos/run119/MNS-R-WT/positions|tss_pm2kb}.bed"
    threads:
        1
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.bed_coordinates} -b {input.bed_feature_to_exclude} -v > {output.bed}"

rule bedtools_intersect_bam_in_feature:
    """
    Created: 2016-03-21
    Modified: 2016-07-21 11h37 - do not take bed from input but from processed. Require some rewriting elsewhere to ensure reproducibility of older tasks. It may be solved by writing a simple rule that create link from input to processed.

    Filter bam for reads falling into specific region types.

    Made to compare the size of insert for reads falling into transcribed genes versus other features.
    expand("out/bedtools/intersect/bam_in_feature/{index}/{exp}/{sample}/{feature}.bam",index=["mm9"],exp=["run125","run126"], sample=["PSK-SC-WT","PSK-SC-KO"], feature=["tss_pm2kb", "microsatellite", "repeatMasker", "segmentalDups", "simpleRepeats","danpos_merge_run113_run119_run124_MNS-R-WT_positions"]),\


    """
    input:
        bam="out/bam/{index}/{exp}/{sample}.bam",
        bed="annotation/processed/feature/{index}/{feature}.bed"
    output:
        bam="out/bedtools/intersect/bam_in_feature/{index, mm[0-9]+}/{exp, [a-zA-Z0-9_-]+}/{sample, [a-zA-Z0-9_-]+}/{feature}.bam"
    threads:
        1
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools intersect -abam {input.bam} -b {input.bed} -wa -u -f 0.5 > {output.bam}"

rule bedtools_intersect_tss_with_0_or_4_PTMs:
    """
    Created: 2016-07-27 10h04 - Take classes with 0 or 4 PTMs and do the heatmaps to see if small structures positions are correlated with PTMs.
    """
    input:
        bed_gene="annotation/input/feature/mm9/knownGene.bed",
        bed_tss="annotation/input/feature/mm9/tss_pm2kb.bed",
        bed_ptm=expand("inp/macs14/{{stage}}-{ptm}_control_{{stage}}-Input/{{stage}}-{ptm}_control_{{stage}}-Input_peaks.bed", ptm=["H4K5ac","H4K5bu","H4K8ac","H4K8bu"])
    output:
        bed_4ptm="out/bedtools/intersect/tss_with_0_or_4_PTMs/{stage, S|R}_4ptm.bed",
        bed_0ptm="out/bedtools/intersect/tss_with_0_or_4_PTMs/{stage, S|R}_0ptm.bed"
    threads:
        1
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        cat {input.bed_gene} > {output.bed_4ptm}_tmp.bed
        cat {input.bed_gene} > {output.bed_0ptm}_tmp.bed

        for bed_ptm in {input.bed_ptm};
        do
            # Finding tss region having one given ptm
            bedtools intersect -a {input.bed_tss} -b $bed_ptm -wa -u > ${{bed_ptm}}_tss.bed
            # Finding genes that have the already looped ptms.
            bedtools intersect -a {output.bed_4ptm}_tmp.bed -b ${{bed_ptm}}_tss.bed -wa -u > {output.bed_4ptm}_tmp2.bed
            cat {output.bed_4ptm}_tmp2.bed > {output.bed_4ptm}_tmp.bed
            # Finding genes that have not the already looped ptms.
            bedtools intersect -a {output.bed_0ptm}_tmp.bed -b ${{bed_ptm}}_tss.bed -wa -v > {output.bed_0ptm}_tmp2.bed
            cat {output.bed_0ptm}_tmp2.bed > {output.bed_0ptm}_tmp.bed
        done

        cat {output.bed_4ptm}_tmp.bed > {output.bed_4ptm}
        cat {output.bed_0ptm}_tmp.bed > {output.bed_0ptm}
        rm -f {output.bed_4ptm}_tmp.bed {output.bed_4ptm}_tmp2.bed {output.bed_0ptm}_tmp.bed {output.bed_0ptm}_tmp2.bed
        """

rule bedtools_intersect_danpos_positions_old:
    """
    Created: 2016-08-30 9h12 - To check this:
    "Les petites structures proviennent majoritairement des positions des nucléosomes de P et R, et minoritairement de SC"

    Control:
    "out/bedtools/shuffle/seed{seed}/{filler}.bed"

    """
    input:
        bed_1="annotation/processed/feature/{index}/danpos_{exp1}_{sample1}_{ppr}.bed",
        bed_2="annotation/processed/feature/{index}/danpos_{exp2}_{sample2}_{ppr}.bed"
    output:
        bed="out/bedtools/intersect/danpos_{ppr, positions|regions|peaks}/{index, mm9|mm10}/{exp1}_{sample1}_VS_{exp2}_{sample2}.bed"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.bed_1} -b {input.bed_2} -f 0.5 -u > {output.bed}"

EXP1=["merge_run113_run125_run126"]*6                                              + ["merge_run113_run119_run124"]*4 + ["merge_run113_run125_run126"]*2
SAMPLE1=["PSK-SC-WT", "PSK-SC-KO"]*3                                               + ["MNS-P-WT", "MNS-P-KO", "MNS-R-WT", "MNS-R-KO", "MNS-SC-WT", "MNS-SC-KO"]
EXP2=["merge_run113_run119_run124"]*4 + ["merge_run113_run125_run126"]*2           + ["merge_run113_run125_run126"]*6
SAMPLE2=["MNS-P-WT", "MNS-P-KO", "MNS-R-WT", "MNS-R-KO", "MNS-SC-WT", "MNS-SC-KO"] + ["PSK-SC-WT","PSK-SC-KO"]*3

rule wc_stats_on_position_intersections_old:
    """
    Created: 2016-08-30 9h12 - To check this:
    "Les petites structures proviennent majoritairement des positions des nucléosomes de P et R, et minoritairement de SC"
    """
    input:
        bed_danpos=expand("annotation/processed/feature/{{index}}/danpos_{exp}_{sample}_{ppr}.bed", zip, ppr=["positions"]*8, exp=["merge_run113_run119_run124"]*4+["merge_run113_run125_run126"]*4, sample=["MNS-P-WT", "MNS-P-KO", "MNS-R-WT", "MNS-R-KO", "MNS-SC-WT", "MNS-SC-KO", "PSK-SC-WT", "PSK-SC-KO"]),
        bed_intersect=expand("out/bedtools/intersect/danpos_{ppr}/{{index}}/{exp1}_{sample1}_VS_{exp2}_{sample2}.bed", zip, ppr=["positions"]*12, exp1=EXP1, sample1=SAMPLE1, exp2=EXP2, sample2=SAMPLE2)
    output:
        wc="result/wc/stats_on_position_intersections/{index}/wc.txt"
    shell:
        """
        wc -l {input.bed_danpos} {input.bed_intersect} > {output.wc}
        """

rule bedtools_intersect_danpos_positions_shuffled:
    """
    Modified: 2016-12-05 15h04 - New pattern to do control with shuffled positions.
    To check this:
    "Les petites structures proviennent majoritairement des positions des nucléosomes de P et R, et minoritairement de SC"

    Control:
    "out/bedtools/shuffle/seed{seed}/danpos/dtriple_v2/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run167/MNase_Spm_WT_band1/danpos.smooth.positions.bed"

    The way seeds work in bedtools shuffle require to have two differents seeds for the samples.
    """
    input:
        bed1="out/bedtools/shuffle/seed{seed1}/{filler}/{exp1}/{sample1}/danpos.smooth.{ppr}.bed",
        bed2="out/bedtools/shuffle/seed{seed2}/{filler}/{exp2}/{sample2}/danpos.smooth.{ppr}.bed"
    output:
        bed="out/bedtools/intersect/danpos_{ppr}_shuffled_seeds_{seed1}_{seed2}/{filler}/{exp1}_{sample1}_VS_{exp2}_{sample2}.bed"
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        bedtools intersect -a {input.bed1} -b {input.bed2} -f 0.5 -u > {output.bed}
        """

rule bedtools_intersect_default:
    """
    Created: 2016-12-16 11h56 - Default usage of bedtools.

    Test:
        "out/bedtools/intersect/default/bed1/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run119_run124/MNS-R-WT/danpos.smooth.positions/bed2/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126/PSK-SC-WT/danpos.smooth.positions.bed"

expand("out/bedtools/genomecov/bed/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run119_run124/{sample}/danpos.smooth.positions.txt", sample=["MNS-P-WT","MNS-P-KO","MNS-R-WT","MNS-R-KO"])$
        expand("out/bedtools/genomecov/bed/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126/{sample}/danpos.smooth.positions.txt", sample=["MNS-SC-WT","MNS-SC-KO","PSK-SC-WT","PSK-SC-KO"])
    """
    input:
        bedtools="opt/bedtools2/bin/bedtools",
        bed1="out/{filler1}.{bedlike}",
        bed2="out/{filler2}.{bedlike}"
    output:
        bed="out/bedtools/intersect/default/bed1/{filler1}/bed2/{filler2}.{bedlike}"
    wildcard_constraints:
        bedlike="bed|bedgraph|gff|vcf"
    shell:
        """
        {input.bedtools} intersect -a {input.bed1} -b {input.bed2} > {output.bed}
        """


rule bedtools_intersect_for_article_indep_rule_with_shuffle:
    """
    Created: 2016-12-14 15h30 - Bedtools shuffle produces strange results with --incl and --seed so I write here a rule that does the whole analysis.
    Bed1 is the bed of the psk and bed2 the one of the mns.

    Test:
        "out/bedtools/intersect/for_article_indep_rule_with_shuffle/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-WT_VS_run113_run119_run124_MNS-P-WT.bed",
        "out/bedtools/intersect/for_article_indep_rule_with_shuffle/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-WT_VS_run113_run119_run124_MNS-R-WT.bed",
        "out/bedtools/intersect/for_article_indep_rule_with_shuffle/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-WT_VS_run113_run125_run126_MNS-SC-WT.bed",
        "out/bedtools/intersect/for_article_indep_rule_with_shuffle/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-KO_VS_run113_run119_run124_MNS-P-KO.bed",
        "out/bedtools/intersect/for_article_indep_rule_with_shuffle/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-KO_VS_run113_run119_run124_MNS-R-KO.bed",
        "out/bedtools/intersect/for_article_indep_rule_with_shuffle/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-KO_VS_run113_run125_run126_MNS-SC-KO.bed"
    """
    input:
        bedtools="opt/bedtools2/bin/bedtools",
        bed1="out/{filler}/{exp1}/{sample1}/danpos.smooth.positions.bed",
        bed2="out/{filler}/{exp2}/{sample2}/danpos.smooth.positions.bed",
        chromInfo="out/mysql/ucsc/mm9/chromInfo_main.tsv"
    output:
        bed_inter="out/bedtools/intersect/for_article_indep_rule_with_shuffle/{filler}/{exp1}_{sample1}_VS_{exp2}_{sample2}.bed",
        bed_inter_shuffle1="out/bedtools/intersect/for_article_indep_rule_with_shuffle/{filler}/{exp1}_{sample1}_VS_{exp2}_{sample2}_shuffle1.bed",
        bed_inter_shuffle2="out/bedtools/intersect/for_article_indep_rule_with_shuffle/{filler}/{exp1}_{sample1}_VS_{exp2}_{sample2}_shuffle2.bed",
        bed_inter_shuffle3="out/bedtools/intersect/for_article_indep_rule_with_shuffle/{filler}/{exp1}_{sample1}_VS_{exp2}_{sample2}_shuffle3.bed",
        bed_inter_shuffle4="out/bedtools/intersect/for_article_indep_rule_with_shuffle/{filler}/{exp1}_{sample1}_VS_{exp2}_{sample2}_shuffle4.bed",
        bed_inter_shuffle5="out/bedtools/intersect/for_article_indep_rule_with_shuffle/{filler}/{exp1}_{sample1}_VS_{exp2}_{sample2}_shuffle5.bed",
        bed_inter_shuffle6="out/bedtools/intersect/for_article_indep_rule_with_shuffle/{filler}/{exp1}_{sample1}_VS_{exp2}_{sample2}_shuffle6.bed",
    params:
        bed2_shuffle1="out/bedtools/intersect/for_article_indep_rule_with_shuffle/{filler}/{exp2}_{sample2}_shuffle1.bed",
        bed2_shuffle2="out/bedtools/intersect/for_article_indep_rule_with_shuffle/{filler}/{exp2}_{sample2}_shuffle2.bed",
        bed2_shuffle3="out/bedtools/intersect/for_article_indep_rule_with_shuffle/{filler}/{exp2}_{sample2}_shuffle3.bed",
        bed2_shuffle4="out/bedtools/intersect/for_article_indep_rule_with_shuffle/{filler}/{exp2}_{sample2}_shuffle4.bed",
        bed2_shuffle5="out/bedtools/intersect/for_article_indep_rule_with_shuffle/{filler}/{exp2}_{sample2}_shuffle5.bed",
        bed2_shuffle6="out/bedtools/intersect/for_article_indep_rule_with_shuffle/{filler}/{exp2}_{sample2}_shuffle6.bed",
    shell:"""
    {input.bedtools} intersect -a {input.bed1} -b {input.bed2} -f 0.5 -u > {output.bed_inter}

    # Shuffle 1
    # Think about that but I think if I shuffle only one of the bed it's enough to produce negative control. I keep small structure bed unchanged so I can shuffle nucleosomes with --noOverlaping.
    {input.bedtools} shuffle \
        -i {input.bed2} \
        -g {input.chromInfo} \
        -noOverlapping > {params.bed2_shuffle1}

    {input.bedtools} intersect -a {input.bed1} -b {params.bed2_shuffle1} -f 0.5 -u > {output.bed_inter_shuffle1}

    # Shuffle 2
    # Think about that but I think if I shuffle only one of the bed it's enough to produce negative control. I keep small structure bed unchanged so I can shuffle nucleosomes with --noOverlaping.
    {input.bedtools} shuffle \
        -i {input.bed2} \
        -g {input.chromInfo} \
        -noOverlapping > {params.bed2_shuffle2}

    {input.bedtools} intersect -a {input.bed1} -b {params.bed2_shuffle2} -f 0.5 -u > {output.bed_inter_shuffle2}

    # Shuffle 3
    # Think about that but I think if I shuffle only one of the bed it's enough to produce negative control. I keep small structure bed unchanged so I can shuffle nucleosomes with --noOverlaping.
    {input.bedtools} shuffle \
        -i {input.bed2} \
        -g {input.chromInfo} \
        -noOverlapping > {params.bed2_shuffle3}

    {input.bedtools} intersect -a {input.bed1} -b {params.bed2_shuffle3} -f 0.5 -u > {output.bed_inter_shuffle3}

    # Shuffle 4
    # Think about that but I think if I shuffle only one of the bed it's enough to produce negative control. I keep small structure bed unchanged so I can shuffle nucleosomes with --noOverlaping.
    {input.bedtools} shuffle \
        -i {input.bed2} \
        -g {input.chromInfo} \
        -noOverlapping > {params.bed2_shuffle4}

    {input.bedtools} intersect -a {input.bed1} -b {params.bed2_shuffle4} -f 0.5 -u > {output.bed_inter_shuffle4}

    # Shuffle 5
    # Think about that but I think if I shuffle only one of the bed it's enough to produce negative control. I keep small structure bed unchanged so I can shuffle nucleosomes with --noOverlaping.
    {input.bedtools} shuffle \
        -i {input.bed2} \
        -g {input.chromInfo} \
        -noOverlapping > {params.bed2_shuffle5}

    {input.bedtools} intersect -a {input.bed1} -b {params.bed2_shuffle5} -f 0.5 -u > {output.bed_inter_shuffle5}

    # Shuffle 6
    # Think about that but I think if I shuffle only one of the bed it's enough to produce negative control. I keep small structure bed unchanged so I can shuffle nucleosomes with --noOverlaping.
    {input.bedtools} shuffle \
        -i {input.bed2} \
        -g {input.chromInfo} \
        -noOverlapping > {params.bed2_shuffle6}

    {input.bedtools} intersect -a {input.bed1} -b {params.bed2_shuffle6} -f 0.5 -u > {output.bed_inter_shuffle6}
    """

rule wc_for_article_noMap:
    input:
        "out/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126/PSK-SC-WT/danpos.smooth.positions.bed",
        "out/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126/PSK-SC-KO/danpos.smooth.positions.bed",
        "out/bedtools/intersect/for_article_indep_rule_with_shuffle/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-WT_VS_run113_run119_run124_MNS-P-WT.bed",
        "out/bedtools/intersect/for_article_indep_rule_with_shuffle/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-WT_VS_run113_run119_run124_MNS-R-WT.bed",
        "out/bedtools/intersect/for_article_indep_rule_with_shuffle/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-WT_VS_run113_run125_run126_MNS-SC-WT.bed",
        "out/bedtools/intersect/for_article_indep_rule_with_shuffle/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-KO_VS_run113_run119_run124_MNS-P-KO.bed",
        "out/bedtools/intersect/for_article_indep_rule_with_shuffle/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-KO_VS_run113_run119_run124_MNS-R-KO.bed",
        "out/bedtools/intersect/for_article_indep_rule_with_shuffle/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-KO_VS_run113_run125_run126_MNS-SC-KO.bed"
    output:
        wc="result/wc/for_article_noMap/wc.txt"
    shell:"""
    wc -l out/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126/PSK-SC-WT/danpos.smooth.positions.bed out/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126/PSK-SC-KO/danpos.smooth.positions.bed out/bedtools/intersect/for_article_indep_rule_with_shuffle/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/* > {output.wc}
    """


rule bedtools_intersect_danpos_positions_shuffled_with_shrink_for_psk:
    """
    Modified: 2016-12-05 15h04 - New pattern to do control with shuffled positions.
    To check this:
    "Les petites structures proviennent majoritairement des positions des nucléosomes de P et R, et minoritairement de SC"

    Control:
    "out/bedtools/shuffle/seed{seed}/danpos/dtriple_v2/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run167/MNase_Spm_WT_band1/danpos.smooth.positions.bed"

    The way seeds work in bedtools shuffle require to have two differents seeds for the samples.
    """
    input:
        bedtools="opt/bedtools2/bin/bedtools",
        bed1="out/bedtools/shuffle/seed{seed1}/awk/shrink_bed/35bp/{filler}/{exp1}/{sample1}/danpos.smooth.{ppr}.bed",
        bed2="out/bedtools/shuffle/seed{seed2}/{filler}/{exp2}/{sample2}/danpos.smooth.{ppr}.bed"
    output:
        bed="out/bedtools/intersect/danpos_{ppr}_shuffled_seeds_{seed1}_{seed2}_with_shrink_for_psk/{filler}/{exp1}_{sample1}_VS_{exp2}_{sample2}.bed"
    shell:"""
    {input.bedtools} intersect -a {input.bed1} -b {input.bed2} -f 0.5 -u > {output.bed}
    """

rule bedtools_intersect_filter_danpos_ppr_on_mapability1:
    """
    Created: 2016-12-08 10h03 - Because I decided to shuffle positions on regions with mapability 1 I think it is more accurate to consider only those positions inside mapability 1 regions.
    "out/r/select_representative_nucleosomes/{filler}/{crit}_min_{npos}.bed"

    Test:
    input:"out/awk/shrink_bed/35bp/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126/PSK-SC-WT/danpos.smooth.positions.bed"
    TODO NEXT: EDIT THIS AND LAUNCH:
    output:"out/bedtools/intersect/filter_danpos_ppr_on_mapability1/awk/shrink_bed/35bp/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126/PSK-SC-WT/danpos.smooth.positions.bed"
    """
    input:
        bedtools="opt/bedtools2/bin/bedtools",
        bed="out/{filler}/{exp}/{sample}/danpos.smooth.{ppr}.bed",
        incl="out/awk/prepare_include_for_shuffle/ucsc/bigWigToBedGraph/rsync/ucsc/gbdb/mm9/bbi/crgMapabilityAlign75mer.bedGraph"
    output:
        bed="out/bedtools/intersect/filter_danpos_ppr_on_mapability1/{filler}/{exp}/{sample}/danpos.smooth.{ppr}.bed"
    shell:"""
    {input.bedtools} intersect -a {input.bed} -b {input.incl} -f 1 -wa > {output.bed}
    """

rule bedtools_intersect_danpos_positions_on_mapability1:
    """
    Created: 2016-12-08 10h03 - Because I decided to shuffle positions on regions with mapability 1 I think it is more accurate to consider only those positions inside mapability 1 regions.
    """
    input:
        bedtools="opt/bedtools2/bin/bedtools",
        bed1="out/bedtools/intersect/filter_danpos_ppr_on_mapability1/{filler}/{exp1}/{sample1}/danpos.smooth.{ppr}.bed",
        bed2="out/bedtools/intersect/filter_danpos_ppr_on_mapability1/{filler}/{exp2}/{sample2}/danpos.smooth.{ppr}.bed"
    output:
        bed="out/bedtools/intersect/danpos_{ppr}_on_mapability1/{filler}/{exp1}_{sample1}_VS_{exp2}_{sample2}.bed"
    shell:"""
    {input.bedtools} intersect -a {input.bed1} -b {input.bed2} -f 0.5 -u > {output.bed}
    """

rule bedtools_intersect_danpos_positions_on_mapability1_with_shrink_for_psk:
    """
    Created: 2016-12-08 10h03 - Because I decided to shuffle positions on regions with mapability 1 I think it is more accurate to consider only those positions inside mapability 1 regions.


    Test:
        output:"out/bedtools/intersect/danpos_positions_on_mapability1_with_shrink_for_psk/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-WT_VS_run113_run119_run124_MNS-P-WT.bed"
    """
    input:
        bedtools="opt/bedtools2/bin/bedtools",
        bed1="out/bedtools/intersect/filter_danpos_ppr_on_mapability1/awk/shrink_bed/35bp/{filler}/{exp1}/{sample1}/danpos.smooth.{ppr}.bed",
        bed2="out/bedtools/intersect/filter_danpos_ppr_on_mapability1/{filler}/{exp2}/{sample2}/danpos.smooth.{ppr}.bed"
    output:
        bed="out/bedtools/intersect/danpos_{ppr}_on_mapability1_with_shrink_for_psk/{filler}/{exp1}_{sample1}_VS_{exp2}_{sample2}.bed"
    shell:"""
    {input.bedtools} intersect -a {input.bed1} -b {input.bed2} -f 0.5 -u > {output.bed}
    """




rule tmp_force_test_mapability1:
    input:
        "out/bedtools/intersect/danpos_positions_on_mapability1/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-WT_VS_run113_run119_run124_MNS-P-WT.bed"


rule wc_stats_on_position_intersections_new_tmp:
    """
    Created: 2016-12-07 11h01 - Stats for shuffled controls.
    """
    input:
        #expand("out/bedtools/intersect/danpos_positions_shuffled_seeds_{seed1}_{seed2}/bedtools/intersect/filter_danpos_ppr_on_mapability1/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-WT_VS_run113_run119_run124_MNS-P-WT.bed", seed1=["1","2","3","4","5"], seed2=["1","2","3","4","5"]),
        #"out/bedtools/intersect/danpos_positions_on_mapability1/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-WT_VS_run113_run119_run124_MNS-P-WT.bed",
        #"out/bedtools/intersect/danpos_positions/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-WT_VS_run113_run119_run124_MNS-P-WT.bed",
        #"out/bedtools/intersect/filter_danpos_ppr_on_mapability1/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126/PSK-SC-WT/danpos.smooth.positions.bed",
        #"out/bedtools/intersect/filter_danpos_ppr_on_mapability1/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run119_run124/MNS-P-WT/danpos.smooth.positions.bed"
        ss_wt="out/bedtools/intersect/filter_danpos_ppr_on_mapability1/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126/PSK-SC-WT/danpos.smooth.positions.bed",
        ss_wt_vs_mns_p="out/bedtools/intersect/danpos_positions_on_mapability1/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-WT_VS_run113_run119_run124_MNS-P-WT.bed",
        ss_wt_vs_mns_p_shuffling=expand("out/bedtools/intersect/danpos_positions_shuffled_seeds_{seed1}_{seed2}/bedtools/intersect/filter_danpos_ppr_on_mapability1/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-WT_VS_run113_run119_run124_MNS-P-WT.bed", zip, seed1=["1","2","3","4","5"], seed2=["5","1","2","3","4"]),
        ss_wt_vs_mns_r="out/bedtools/intersect/danpos_positions_on_mapability1/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-WT_VS_run113_run119_run124_MNS-R-WT.bed",
        ss_wt_vs_mns_r_shuffling=expand("out/bedtools/intersect/danpos_positions_shuffled_seeds_{seed1}_{seed2}/bedtools/intersect/filter_danpos_ppr_on_mapability1/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-WT_VS_run113_run119_run124_MNS-R-WT.bed", zip, seed1=["1","2","3","4","5"], seed2=["5","1","2","3","4"]),
        ss_wt_vs_mns_sc="out/bedtools/intersect/danpos_positions_on_mapability1/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-WT_VS_run113_run125_run126_MNS-SC-WT.bed",
        ss_wt_vs_mns_sc_shuffling=expand("out/bedtools/intersect/danpos_positions_shuffled_seeds_{seed1}_{seed2}/bedtools/intersect/filter_danpos_ppr_on_mapability1/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-WT_VS_run113_run125_run126_MNS-SC-WT.bed", zip, seed1=["1","2","3","4","5"], seed2=["5","1","2","3","4"]),
        ss_ko="out/bedtools/intersect/filter_danpos_ppr_on_mapability1/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126/PSK-SC-KO/danpos.smooth.positions.bed",
        ss_ko_vs_mns_p="out/bedtools/intersect/danpos_positions_on_mapability1/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-KO_VS_run113_run119_run124_MNS-P-KO.bed",
        ss_ko_vs_mns_p_shuffling=expand("out/bedtools/intersect/danpos_positions_shuffled_seeds_{seed1}_{seed2}/bedtools/intersect/filter_danpos_ppr_on_mapability1/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-KO_VS_run113_run119_run124_MNS-P-KO.bed", zip, seed1=["1","2","3","4","5"], seed2=["5","1","2","3","4"]),
        ss_ko_vs_mns_r="out/bedtools/intersect/danpos_positions_on_mapability1/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-KO_VS_run113_run119_run124_MNS-R-KO.bed",
        ss_ko_vs_mns_r_shuffling=expand("out/bedtools/intersect/danpos_positions_shuffled_seeds_{seed1}_{seed2}/bedtools/intersect/filter_danpos_ppr_on_mapability1/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-KO_VS_run113_run119_run124_MNS-R-KO.bed", zip, seed1=["1","2","3","4","5"], seed2=["5","1","2","3","4"]),
        ss_ko_vs_mns_sc="out/bedtools/intersect/danpos_positions_on_mapability1/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-KO_VS_run113_run125_run126_MNS-SC-KO.bed",
        ss_ko_vs_mns_sc_shuffling=expand("out/bedtools/intersect/danpos_positions_shuffled_seeds_{seed1}_{seed2}/bedtools/intersect/filter_danpos_ppr_on_mapability1/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-KO_VS_run113_run125_run126_MNS-SC-KO.bed", zip, seed1=["1","2","3","4","5"], seed2=["5","1","2","3","4"]),

    output:
        wc="result/wc/stats_on_position_intersections_new_tmp/{index}/wc.txt"
    shell:"""
    wc -l {input} > {output.wc}
    """

rule wc_stats_on_position_intersections_with_shrink_for_psk:
    """
    Created: 2016-12-13 15h45 - Integrating the concept of reducing small structures position size.
    """
    input:
        ss_wt="out/bedtools/intersect/filter_danpos_ppr_on_mapability1/awk/shrink_bed/35bp/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126/PSK-SC-WT/danpos.smooth.positions.bed",
        ss_wt_vs_mns_p="out/bedtools/intersect/danpos_positions_on_mapability1_with_shrink_for_psk/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-WT_VS_run113_run119_run124_MNS-P-WT.bed",
        ss_wt_vs_mns_p_shuffling=expand("out/bedtools/intersect/danpos_positions_shuffled_seeds_{seed1}_{seed2}_with_shrink_for_psk/bedtools/intersect/filter_danpos_ppr_on_mapability1/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-WT_VS_run113_run119_run124_MNS-P-WT.bed", zip, seed1=["1","2","3","4","5"], seed2=["5","1","2","3","4"]),
        ss_wt_vs_mns_r="out/bedtools/intersect/danpos_positions_on_mapability1_with_shrink_for_psk/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-WT_VS_run113_run119_run124_MNS-R-WT.bed",
        ss_wt_vs_mns_r_shuffling=expand("out/bedtools/intersect/danpos_positions_shuffled_seeds_{seed1}_{seed2}_with_shrink_for_psk/bedtools/intersect/filter_danpos_ppr_on_mapability1/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-WT_VS_run113_run119_run124_MNS-R-WT.bed", zip, seed1=["1","2","3","4","5"], seed2=["5","1","2","3","4"]),
        ss_wt_vs_mns_sc="out/bedtools/intersect/danpos_positions_on_mapability1_with_shrink_for_psk/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-WT_VS_run113_run125_run126_MNS-SC-WT.bed",
        ss_wt_vs_mns_sc_shuffling=expand("out/bedtools/intersect/danpos_positions_shuffled_seeds_{seed1}_{seed2}_with_shrink_for_psk/bedtools/intersect/filter_danpos_ppr_on_mapability1/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126_PSK-SC-WT_VS_run113_run125_run126_MNS-SC-WT.bed", zip, seed1=["1","2","3","4","5"], seed2=["5","1","2","3","4"]),
    output:
        wc="result/wc/stats_on_position_intersections_with_shrink_for_psk/mm9/wc.txt"
    shell:"""
    wc -l {input} > {output.wc}
    """

rule bedtools_intersect_tss_with_brdt_and_h4k5ac:
    """
    Created:
        2017-01-24 16h08
    Test:
        "out/bedtools/intersect/tss_with_brdt_and_h4k5ac.bed"
    """
    input:
        bedtools="opt/bedtools2/bin/bedtools",
        bed_brdt="out/gunzip/to-stdout/tar/xvf_data_brdt/GSM984200_4_EM2_R_BRDT-Paired.F3.no_4-W200-G600-islands-summary-FDR001.bed",
        bed_h4k5ac="out/awk/extract_xls_coordinates_to_bed3/macs2/callpeak_broad_BAMPE_mm_SPMRTrue/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/sickle/pe_-t_sanger_-q_30/merge_lanes/run140_run141/H4K5ac-Nut-WT_over_Input-Nut-WT_peaks.bed",
        bed_tss="out/sort/unique_coordinates_bed/inp/annotation/input/feature/mm9/tss_pm2kb.bed"
    output:
        bed="out/bedtools/intersect/mm9_tss_with_brdt_and_h4k5ac.bed"
    #missing zlib in conda env.
    #conda:
    #    "../envs/bedtools.yaml"
    shell:
        """
        {input.bedtools} intersect -a {input.bed_tss} -b {input.bed_brdt} -u -wa -f 0.5 | \
            {input.bedtools} intersect -a stdin -b {input.bed_h4k5ac} -u -wa -f 0.5 > {output.bed}
        """

rule bedtools_intersect_tss_with_brd4_and_h4k5ac:
    """
    Created: 2017-01-30 13h25

    Test:
        "out/bedtools/intersect/tss_with_brd4_and_h4k5ac.bed"
    """
    input:
        bedtools="opt/bedtools2/bin/bedtools",
        bed_brd4="out/awk/extract_xls_coordinates_to_bed3/macs2/callpeak_broad_BAM_mm_SPMRTrue/samtools/sam_to_bam/bowtie2/se/mm9/GSE56526/SRR1596612_over_SRR1596620_peaks.bed",
        bed_h4k5ac="out/awk/extract_xls_coordinates_to_bed3/macs2/callpeak_broad_BAMPE_mm_SPMRTrue/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/sickle/pe_-t_sanger_-q_30/merge_lanes/run140_run141/H4K5ac-Nut-WT_over_Input-Nut-WT_peaks.bed",
        #bed_h4k5ac="out/awk/extract_xls_coordinates_to_bed3/macs2/callpeak_broad_BAMPE_mm_SPMRTrue/bam/mm9/merge_run140_run141/H4K5ac-Nut-WT_over_Input-Nut-WT_peaks.bed",
        bed_tss="out/sort/unique_coordinates_bed/inp/annotation/input/feature/mm9/tss_pm2kb.bed"
        #bed_tss="out/sort/unique_coordinates_bed/input/annotation/feature/mm9/tss_pm2kb.bed"
    output:
        bed="out/bedtools/intersect/mm9_tss_with_brd4_and_h4k5ac.bed"
    #missing zlib in conda env.
    #conda:
    #    "../envs/bedtools.yaml"
    shell:
        """
        {input.bedtools} intersect -a {input.bed_tss} -b {input.bed_brd4} -u -wa -f 0.3 | \
            {input.bedtools} intersect -a stdin -b {input.bed_h4k5ac} -u -wa -f 0.5 > {output.bed}
        """

rule bedtools_intersect_tss_without_brdt_brd4_and_with_h4k5ac:
    """
    Created: 2017-01-30 13h25
    """
    input:
        bedtools="opt/bedtools2/bin/bedtools",
        bed_brdt="out/gunzip/tar/xvf/data_brdt/GSM984200_4_EM2_R_BRDT-Paired.F3.no_4-W200-G600-islands-summary-FDR001.bed",
        bed_brd4="out/awk/extract_xls_coordinates_to_bed3/macs2/callpeak_broad_BAM_mm_SPMRTrue/samtools/sam_to_bam/bowtie2/se/mm9/GSE56526/SRR1596612_over_SRR1596620_peaks.bed",
        bed_h4k5ac="out/awk/extract_xls_coordinates_to_bed3/macs2/callpeak_broad_BAMPE_mm_SPMRTrue/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/sickle/pe_-t_sanger_-q_30/merge_lanes/run140_run141/H4K5ac-Nut-WT_over_Input-Nut-WT_peaks.bed",
        bed_tss="out/sort/unique_coordinates_bed/inp/annotation/input/feature/mm9/tss_pm2kb.bed"
        #bed_h4k5ac="out/awk/extract_xls_coordinates_to_bed3/macs2/callpeak_broad_BAMPE_mm_SPMRTrue/bam/mm9/merge_run140_run141/H4K5ac-Nut-WT_over_Input-Nut-WT_peaks.bed",
        #bed_tss="out/sort/unique_coordinates_bed/input/annotation/feature/mm9/tss_pm2kb.bed"
    output:
        bed="out/bedtools/intersect/mm9_tss_without_brdt_brd4_and_with_h4k5ac.bed"
    #missing zlib in conda env.
    #conda:
    #    "../envs/bedtools.yaml"
    shell:
        """
        {input.bedtools} intersect -a {input.bed_tss} -b {input.bed_brd4} -v -wa -f 0.1 | \
            {input.bedtools} intersect -a stdin -b {input.bed_brdt} -v -wa -f 0.1 | \
            {input.bedtools} intersect -a stdin -b {input.bed_h4k5ac} -u -wa -f 0.5 > {output.bed}
        """

rule bedtools_intersect_tss_with_h4k5ac:
    """
    Created:
        2017-03-09 15:34:22
    Aim:
        Get a reference for coverage for tss with H4K5ac to compare with classes with Brdt/Brd4.
    Note:
        In these rules I used to ask a minimum overlap of 0.1 fraction but I think I may get better result if I increase the fraction in common. So at 2017-03-09 15:42:40 I modified all inclusions fraction requirement to 0.6 and kept exlclusion to 0.1.
    """
    input:
        bedtools="opt/bedtools2/bin/bedtools",
        bed_h4k5ac="out/awk/extract_xls_coordinates_to_bed3/macs2/callpeak_broad_BAMPE_mm_SPMRTrue/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/sickle/pe_-t_sanger_-q_30/merge_lanes/run140_run141/H4K5ac-Nut-WT_over_Input-Nut-WT_peaks.bed",
        bed_tss="out/sort/unique_coordinates_bed/inp/annotation/input/feature/mm9/tss_pm2kb.bed"
        #bed_h4k5ac="out/awk/extract_xls_coordinates_to_bed3/macs2/callpeak_broad_BAMPE_mm_SPMRTrue/bam/mm9/merge_run140_run141/H4K5ac-Nut-WT_over_Input-Nut-WT_peaks.bed",
        #bed_tss="out/sort/unique_coordinates_bed/input/annotation/feature/mm9/tss_pm2kb.bed"
    output:
        bed="out/bedtools/intersect/mm9_tss_with_h4k5ac.bed"
    shell:
        """
        {input.bedtools} intersect -a {input.bed_tss} -b {input.bed_h4k5ac} -u -wa -f 0.5 > {output.bed}
        """

rule bedtools_intersect_tss_nut_updown_with_brdt:
    """
    Created:
        2017-03-23 11:54:45
    Aim:
        Look at fraction of genes up/down-regulated that have Brdt.
    Note:
        ("out/sort/unique_coordinates_bed/awk/extract_gfold_signif_gene_threshold1/star/pe_mm9/sickle/pe_-t_sanger_-q_20/fastq/nut_rnaseq/Nut-R-WT_VS_Nut-R-KO_{posneg}.bed", posneg=["pos","neg"])
    Test:
        out/bedtools/intersect/mm9_tss_nut_updown_with_brdt_pos.bed
    """
    input:
        bedtools="opt/bedtools2/bin/bedtools",
        bed_tss_class="out/sort/unique_coordinates_bed/awk/extract_gfold_signif_gene_threshold1/star/pe_mm9/sickle/pe_-t_sanger_-q_20/fastq/nut_rnaseq/Nut-R-WT_VS_Nut-R-KO_{pos_neg}.bed",
        bed_brdt="out/gunzip/tar/xvf/data_brdt/GSM984200_4_EM2_R_BRDT-Paired.F3.no_4-W200-G600-islands-summary-FDR001.bed"
    output:
        bed="out/bedtools/intersect/mm9_tss_nut_updown_with_brdt_{pos_neg}.bed"
    shell:
        """
        {input.bedtools} intersect -a {input.bed_tss_class} -b {input.bed_brdt} -u -wa -f 0.5 > {output.bed}
        """

rule bedtools_intersect_tss_nut_updown_without_brdt:
    """
    Created:
        2017-03-23 11:54:45
    Aim:
        Look at fraction of genes up/down-regulated that have Brdt.
    Note:
        ("out/sort/unique_coordinates_bed/awk/extract_gfold_signif_gene_threshold1/star/pe_mm9/sickle/pe_-t_sanger_-q_20/fastq/nut_rnaseq/Nut-R-WT_VS_Nut-R-KO_{posneg}.bed", posneg=["pos","neg"])
    Test:
        out/bedtools/intersect/mm9_tss_nut_updown_without_brdt_neg.bed
    """
    input:
        bedtools="opt/bedtools2/bin/bedtools",
        bed_tss_class="out/sort/unique_coordinates_bed/awk/extract_gfold_signif_gene_threshold1/star/pe_mm9/sickle/pe_-t_sanger_-q_20/fastq/nut_rnaseq/Nut-R-WT_VS_Nut-R-KO_{pos_neg}.bed",
        bed_brdt="out/gunzip/tar/xvf/data_brdt/GSM984200_4_EM2_R_BRDT-Paired.F3.no_4-W200-G600-islands-summary-FDR001.bed"
    output:
        bed="out/bedtools/intersect/mm9_tss_nut_updown_without_brdt_{pos_neg}.bed"
    shell:
        """
        {input.bedtools} intersect -a {input.bed_tss_class} -b {input.bed_brdt} -v -wa -f 0.1 > {output.bed}
        """

rule bedtools_intersect_complementary_features:
    """
    Created:
        2017-04-12 14:29:34
    Aim:
        I have a bed file with regions captured and I want to define the bed file of regions not captured (for G38) in order to use as --blacklist argument in deepTools.(-region argument does not take bed as input).
    Test:
    Note:
        Probably useless to develop it because bedtools complement does what I want to do.
    """
    input:
        bed_coordinates="out/awk/chromInfo_to_bed3/gunzip/rsync/ucsc/goldenPath/hg38/database/chromInfo.bed",
        bed_feature_to_exclude="out/crossmap/bed_hg19_to_hg38/input/bed/salva/Regions_capture_SE_Starr-seq_Alex_HG19_Merged_cellLine_specific_SE.bed",
        bedtools="opt/bedtools2/bin/bedtools"
    output:
        bed="out/bedtools/intersect_complementary_features/crossmap/bed_hg19_to_hg38/input/bed/salva/Regions_capture_SE_Starr-seq_Alex_HG19_Merged_cellLine_specific_SE.bed"
    threads:
        1
    shell:
        """
        {input.bedtools} intersect -a {input.bed_coordinates} -b {input.bed_feature_to_exclude} -v > {output.bed}
        """

rule bedtools_intersect_extract_ss_in_nuc_test1:
    """
    Created:
        2017-02-06 14:42:44
    Aim:
        Extract small structures positions that fall into nucleosomes positions from round spermatid.
    Test:
        out/bedttols/intersect_extract_ss_in_nuc/test1.bed
    """
    input:
        bedtools="opt/bedtools2/bin/bedtools",
        bed_ss="out/awk/shrink_bed/35bp/awk/extract_danpos_pos_to_bed_maxfuz35_minsmt200/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126/PSK-SC-WT/danpos.smooth.positions.bed",
        bed_nuc="out/awk/extract_danpos_pos_to_bed_maxfuz40_minsmt200/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run119_run124/MNS-R-WT/danpos.smooth.positions.bed"
    output:
        bed_ss="out/bedtools/intersect_extract_ss_in_nuc/test1.bed"
    shell:
        """
        {input.bedtools} intersect -a {input.bed_ss} -b {input.bed_nuc} -wa -u -f 0.9 > {output.bed_ss}
        """

rule bedtools_intersect_extract_ss_in_nuc_test2:
    """
    Created:
        2017-02-07 15:02:25
    Aim:
        Extract small structures positions that fall into nucleosomes positions from round spermatid. This one return the original danpos positions for small structure (140bp) instead of the shrinked one to 70bp. I need to do some controls to see if I need to modify the step size in DANPOS caller.
    Test:
        out/bedtools/intersect_extract_ss_in_nuc/test2_ss.bed
    """
    input:
        bedtools="opt/bedtools2/bin/bedtools",
        bed_ss="out/awk/extract_danpos_pos_to_bed_maxfuz35_minsmt200/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126/PSK-SC-WT/danpos.smooth.positions.bed",
        bed_nuc="out/awk/extract_danpos_pos_to_bed_maxfuz40_minsmt200/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run119_run124/MNS-R-WT/danpos.smooth.positions.bed"
    output:
        bed_ss="out/bedtools/intersect_extract_ss_in_nuc/test2_ss.bed",
        bed_nuc="out/bedtools/intersect_extract_ss_in_nuc/test2_nuc.bed"
    shell:
        """
        {input.bedtools} intersect -a {input.bed_ss} -b {input.bed_nuc} -wa -u -f 0.5 > {output.bed_ss}
        {input.bedtools} intersect -a {input.bed_nuc} -b {input.bed_ss} -wa -u -f 0.5 > {output.bed_nuc}
        """

rule bedtools_intersect_exclude_from_blacklist_mm9:
    """
    Created:
        2017-10-27 16:11:21
    Aim:

    Test:
    """
    input:
        bedtools="opt/bedtools2/bin/bedtools",
        bed="out/{filler}.bed",
        bed_blacklist="out/gunzip/to-stdout/wget/mitra_stanford_ed/kundaje/akundaje/release/blacklists/mm9-mouse/mm9-blacklist.bed"
    output:
        bed="out/bedtools/intersect_exclude_from_blacklist_mm9/{filler}.bed"
    shell:
        """
        {input.bedtools} intersect -a {input.bed} -b {input.bed_blacklist} -v > {output.bed}
        """

rule bedtools_intersect_exclude_from_blacklist_mm10:
    """
    Created:
        2017-08-03 10:33:47
    Aim:

    Test:
        out/bedtools/intersect_exclude_from_blacklist_mm10/awk/extract_gfold_signif_peaks_threshold-0/gfold/diff/gfold/count_bed-mm10-R-H4K5ac-bu-peaks/samtools/sort/samtools/view_bSh/bowtie2/se_mm10/sickle/se_-t_sanger_-q_20/sra-tools/fastqdump_se/wget/ftp_trace_ncbi_sra/SRR3126243_VS_SRR3126244_pos.bed
    """
    input:
        bedtools="opt/bedtools2/bin/bedtools",
        bed_mm10="out/{filler}.bed",
        bed_blacklist="out/gunzip/to-stdout/wget/mitra_stanford_ed/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed"
    output:
        bed="out/bedtools/intersect_exclude_from_blacklist_mm10/{filler}.bed"
    shell:
        """
        {input.bedtools} intersect -a {input.bed_mm10} -b {input.bed_blacklist} -v > {output.bed}
        """

rule bedtools_intersect_exclude_from_gene_tss_mm10:
    """
    Created:
        2017-08-03 10:33:47
    Aim:
        We want to see regions with H4K5/8-ac/bu outside TSS.
    Test:

        out/bedtools/intersect_exclude_from_gene_tss_mm10/awk/extract_gfold_signif_peaks_threshold-0/gfold/diff/gfold/count_bed-mm10-R-H4K5ac-bu-peaks/samtools/sort/samtools/view_bSh/bowtie2/se_mm10/sickle/se_-t_sanger_-q_20/sra-tools/fastqdump_se/wget/ftp_trace_ncbi_sra/SRR3126243_VS_SRR3126244_pos.bed
    """
    input:
        bedtools="opt/bedtools2/bin/bedtools",
        bed_mm10="out/{filler}.bed",
        bed_5p_coord="out/gtftk/5p_3p_coord/sed/add_chr/gunzip/to-stdout/wget/ftp_ensembl/pub/release-87/gtf/mus_musculus/Mus_musculus.GRCm38.87.bed"
    output:
        bed="out/bedtools/intersect_exclude_from_gene_tss_mm10/{filler}.bed"
    shell:
        """
        {input.bedtools} intersect -a {input.bed_mm10} -b {input.bed_5p_coord} -v > {output.bed}
        """

rule bedtools_intersect_exclude_from_gene_tss_5kb_mm10:
    """
    Created:
        2017-08-03 13:50:15
    Aim:
        We want to see regions with H4K5/8-ac/bu outside TSS.
    Test:

        out/bedtools/intersect_exclude_from_gene_tss_mm10/awk/extract_gfold_signif_peaks_threshold-0/gfold/diff/gfold/count_bed-mm10-R-H4K5ac-bu-peaks/samtools/sort/samtools/view_bSh/bowtie2/se_mm10/sickle/se_-t_sanger_-q_20/sra-tools/fastqdump_se/wget/ftp_trace_ncbi_sra/SRR3126243_VS_SRR3126244_pos.bed
    """
    input:
        bedtools="opt/bedtools2/bin/bedtools",
        bed_mm10="out/{filler}.bed",
        bed_5p_coord="out/bedtools/slop_b-5000/gtftk/5p_3p_coord/sed/add_chr/gunzip/to-stdout/wget/ftp_ensembl/pub/release-87/gtf/mus_musculus/Mus_musculus.GRCm38.87.bed"
    output:
        bed="out/bedtools/intersect_exclude_from_gene_tss_5kb_mm10/{filler}.bed"
    shell:
        """
        {input.bedtools} intersect -a {input.bed_mm10} -b {input.bed_5p_coord} -v > {output.bed}
        """

rule bedtools_intersect_exclude_from_macs2_peaks_in_input_hg38:
    """
    Created:
        2018-02-02 01:39:53
    Aim:
        Removing peaks from blacklisted regions.
    Test:
        out/bedtools/intersect_exclude_from_macs2_peaks_in_input_hg38/bedtools/merge/sort/coordinates_bed/cat/hg38-macs2-peaks-H3K27ac-thymus.bed
    """
    input:
        bedtools="opt/bedtools2/bin/bedtools",
        bed_a="out/{filler}.bed",
        bed_b="out/macs2/callpeak_broad_format-BAM_gsize-hs_no_control/ln/alias/experiments/hg38_H3K27ac_thymus/input_peaks.broadPeak"
    output:
        bed="out/bedtools/intersect_exclude_from_macs2_peaks_in_input_hg38/{filler}.bed"
    shell:
        """
        {input.bedtools} intersect -a {input.bed_a} -b {input.bed_b} -v > {output.bed}
        """

rule bedtools_intersect_thymus_stage_exclusive_peaks_hg38:
    """
    Created:
        2018-02-13 09:49:05
    Aim:
        Creating thymus stage exclusive peaks from macs2 peaks.
    Test:
        out/bedtools/intersect_exclude_from_macs2_peaks_in_input_hg38/bedtools/intersect_thymus_stage_exclusive_peaks_hg38/CD34_134.bed
        "Th134_CD34", "Th89_EC", "Th101_LC", "Th125_SP4", "Th125_SP8"
    """
    input:
        bedtools="opt/bedtools2/bin/bedtools",
        bed_cd34_134="out/macs2/callpeak_broad_format-BAM_gsize-hs/ln/alias/experiments/{assembly}_H3K27ac_thymus/Th134_CD34_over_input_peaks.broadPeak",
        bed_ec_89   ="out/macs2/callpeak_broad_format-BAM_gsize-hs/ln/alias/experiments/{assembly}_H3K27ac_thymus/Th89_EC_over_input_peaks.broadPeak",
        bed_ec_125  ="out/macs2/callpeak_broad_format-BAM_gsize-hs/ln/alias/experiments/{assembly}_H3K27ac_thymus/Th125_EC_over_input_peaks.broadPeak",
        bed_lc_91   ="out/macs2/callpeak_broad_format-BAM_gsize-hs/ln/alias/experiments/{assembly}_H3K27ac_thymus/Th91_LC_over_input_peaks.broadPeak",
        bed_lc_101  ="out/macs2/callpeak_broad_format-BAM_gsize-hs/ln/alias/experiments/{assembly}_H3K27ac_thymus/Th101_LC_over_input_peaks.broadPeak",
        bed_lc_118  ="out/macs2/callpeak_broad_format-BAM_gsize-hs/ln/alias/experiments/{assembly}_H3K27ac_thymus/Th118_LC_over_input_peaks.broadPeak",
        bed_lc_125  ="out/macs2/callpeak_broad_format-BAM_gsize-hs/ln/alias/experiments/{assembly}_H3K27ac_thymus/Th125_LC_over_input_peaks.broadPeak",
        bed_sp4_91  ="out/macs2/callpeak_broad_format-BAM_gsize-hs/ln/alias/experiments/{assembly}_H3K27ac_thymus/Th91_SP4_over_input_peaks.broadPeak",
        bed_sp4_125 ="out/macs2/callpeak_broad_format-BAM_gsize-hs/ln/alias/experiments/{assembly}_H3K27ac_thymus/Th125_SP4_over_input_peaks.broadPeak",
        bed_sp8_125 ="out/macs2/callpeak_broad_format-BAM_gsize-hs/ln/alias/experiments/{assembly}_H3K27ac_thymus/Th125_SP8_over_input_peaks.broadPeak"
        ### ADD Peaks from each stage. Need to select which peaks we get from each stage.
    output:
        bed_cd34_134="out/bedtools/intersect_thymus_stage_exclusive_peaks_{assembly}/CD34_134.bed",
        bed_ec_89   ="out/bedtools/intersect_thymus_stage_exclusive_peaks_{assembly}/EC_89.bed",
        #bed_ec_125  ="out/bedtools/intersect_thymus_stage_exclusive_peaks_{assembly}/EC_125.bed",
        #bed_lc_91   ="out/bedtools/intersect_thymus_stage_exclusive_peaks_{assembly}/LC_91.bed",
        bed_lc_101  ="out/bedtools/intersect_thymus_stage_exclusive_peaks_{assembly}/LC_101.bed",
        #bed_lc_118  ="out/bedtools/intersect_thymus_stage_exclusive_peaks_{assembly}/LC_118.bed",
        #bed_lc_125  ="out/bedtools/intersect_thymus_stage_exclusive_peaks_{assembly}/LC_125.bed",
        #bed_sp4_91  ="out/bedtools/intersect_thymus_stage_exclusive_peaks_{assembly}/SP4_91.bed"
        bed_sp4_125 ="out/bedtools/intersect_thymus_stage_exclusive_peaks_{assembly}/SP4_125.bed",
        bed_sp8_125 ="out/bedtools/intersect_thymus_stage_exclusive_peaks_{assembly}/SP8_125.bed",
    shell:
        """
        {input.bedtools} intersect -a {input.bed_cd34_134} -b {input.bed_ec_89} {input.bed_lc_101} {input.bed_sp4_125} {input.bed_sp8_125} -v > {output.bed_cd34_134}
        # $ wc -l out/bedtools/intersect_exclude_from_macs2_peaks_in_input_hg38/bedtools/intersect_thymus_stage_exclusive_peaks_hg38/CD34_134.bed          # 28254 out/bedtools/intersect_exclude_from_macs2_peaks_in_input_hg38/bedtools/intersect_thymus_stage_exclusive_peaks_hg38/CD34_134.bed
        #{input.bedtools} intersect -a {input.bed_cd34_134} -b {input.bed_ec_89} -v |\
        #{input.bedtools} intersect -a stdin -b {input.bed_lc_101} -v |\
        #{input.bedtools} intersect -a stdin -b {input.bed_sp4_125} -v |\
        #{input.bedtools} intersect -a stdin -b {input.bed_sp8_125} -v > {output.bed_cd34_134}

        {input.bedtools} intersect -a {input.bed_ec_89} -b {input.bed_cd34_134} {input.bed_lc_101} {input.bed_sp4_125} {input.bed_sp8_125} -v > {output.bed_ec_89}
        {input.bedtools} intersect -a {input.bed_lc_101} -b {input.bed_cd34_134} {input.bed_ec_89} {input.bed_sp4_125} {input.bed_sp8_125} -v > {output.bed_lc_101}
        {input.bedtools} intersect -a {input.bed_sp4_125} -b {input.bed_cd34_134} {input.bed_ec_89} {input.bed_lc_101} {input.bed_sp8_125} -v > {output.bed_sp4_125}
        {input.bedtools} intersect -a {input.bed_sp8_125} -b {input.bed_cd34_134} {input.bed_ec_89} {input.bed_lc_101} {input.bed_sp4_125} -v > {output.bed_sp8_125}
        """

rule bedtools_intersect_thymus_stage_exclusive_peaks_hg38_v3:
    """
    Created:
        2018-02-25 12:56:26
    Aim:
        Creating thymus stage exclusive peaks from macs2 peaks.
        V3 because I did not take the samples Salva wanted in the previous intersection.
    Test:
        out/bedtools/intersect_exclude_from_macs2_peaks_in_input_hg38/bedtools/intersect_thymus_stage_exclusive_peaks_hg38/CD34_134.bed
        "Th134_CD34", "Th89_EC", "Th101_LC", "Th125_SP4", "Th125_SP8"
    """
    input:
        bedtools="opt/bedtools2/bin/bedtools",
        bed_cd34_134="out/macs2/callpeak_broad_format-BAM_gsize-hs/ln/alias/experiments/{assembly}_H3K27ac_thymus/Th134_CD34_over_input_peaks.broadPeak",
        bed_ec_89   ="out/macs2/callpeak_broad_format-BAM_gsize-hs/ln/alias/experiments/{assembly}_H3K27ac_thymus/Th89_EC_over_input_peaks.broadPeak",
        bed_ec_125  ="out/macs2/callpeak_broad_format-BAM_gsize-hs/ln/alias/experiments/{assembly}_H3K27ac_thymus/Th125_EC_over_input_peaks.broadPeak",
        bed_lc_91   ="out/macs2/callpeak_broad_format-BAM_gsize-hs/ln/alias/experiments/{assembly}_H3K27ac_thymus/Th91_LC_over_input_peaks.broadPeak",
        bed_lc_101  ="out/macs2/callpeak_broad_format-BAM_gsize-hs/ln/alias/experiments/{assembly}_H3K27ac_thymus/Th101_LC_over_input_peaks.broadPeak",
        bed_lc_118  ="out/macs2/callpeak_broad_format-BAM_gsize-hs/ln/alias/experiments/{assembly}_H3K27ac_thymus/Th118_LC_over_input_peaks.broadPeak",
        bed_lc_125  ="out/macs2/callpeak_broad_format-BAM_gsize-hs/ln/alias/experiments/{assembly}_H3K27ac_thymus/Th125_LC_over_input_peaks.broadPeak",
        bed_sp4_91  ="out/macs2/callpeak_broad_format-BAM_gsize-hs/ln/alias/experiments/{assembly}_H3K27ac_thymus/Th91_SP4_over_input_peaks.broadPeak",
        bed_sp4_125 ="out/macs2/callpeak_broad_format-BAM_gsize-hs/ln/alias/experiments/{assembly}_H3K27ac_thymus/Th125_SP4_over_input_peaks.broadPeak",
        bed_sp8_125 ="out/macs2/callpeak_broad_format-BAM_gsize-hs/ln/alias/experiments/{assembly}_H3K27ac_thymus/Th125_SP8_over_input_peaks.broadPeak"
        ### ADD Peaks from each stage. Need to select which peaks we get from each stage.
    output:
        bed_cd34_134="out/bedtools/intersect_thymus_stage_exclusive_peaks_{assembly}_v3/CD34_134.bed",
        bed_ec_89   ="out/bedtools/intersect_thymus_stage_exclusive_peaks_{assembly}_v3/EC_89.bed",
        #bed_ec_125  ="out/bedtools/intersect_thymus_stage_exclusive_peaks_{assembly}_v3/EC_125.bed",
        bed_lc_91   ="out/bedtools/intersect_thymus_stage_exclusive_peaks_{assembly}_v3/LC_91.bed",
        #bed_lc_101  ="out/bedtools/intersect_thymus_stage_exclusive_peaks_{assembly}_v3/LC_101.bed",
        #bed_lc_118  ="out/bedtools/intersect_thymus_stage_exclusive_peaks_{assembly}_v3/LC_118.bed",
        #bed_lc_125  ="out/bedtools/intersect_thymus_stage_exclusive_peaks_{assembly}_v3/LC_125.bed",
        bed_sp4_91  ="out/bedtools/intersect_thymus_stage_exclusive_peaks_{assembly}_v3/SP4_91.bed",
        #bed_sp4_125 ="out/bedtools/intersect_thymus_stage_exclusive_peaks_{assembly}_v3/SP4_125.bed",
        bed_sp8_125 ="out/bedtools/intersect_thymus_stage_exclusive_peaks_{assembly}_v3/SP8_125.bed",
    shell:
        """
        {input.bedtools} intersect -a {input.bed_cd34_134} -b {input.bed_ec_89} {input.bed_lc_91} {input.bed_sp4_91} {input.bed_sp8_125} -v > {output.bed_cd34_134}
        {input.bedtools} intersect -a {input.bed_ec_89} -b {input.bed_cd34_134} {input.bed_lc_91} {input.bed_sp4_91} {input.bed_sp8_125} -v > {output.bed_ec_89}
        {input.bedtools} intersect -a {input.bed_lc_91} -b {input.bed_cd34_134} {input.bed_ec_89} {input.bed_sp4_91} {input.bed_sp8_125} -v > {output.bed_lc_91}
        {input.bedtools} intersect -a {input.bed_sp4_91} -b {input.bed_cd34_134} {input.bed_ec_89} {input.bed_lc_91} {input.bed_sp8_125} -v > {output.bed_sp4_91}
        {input.bedtools} intersect -a {input.bed_sp8_125} -b {input.bed_cd34_134} {input.bed_ec_89} {input.bed_lc_91} {input.bed_sp4_91} -v > {output.bed_sp8_125}
        """

rule bedtools_intersect_thymus_stage_exclusive_peaks_hg38_merged:
    """
    Created:
        2018-03-16 19:26:12
    Aim:
        Creating thymus stage exclusive peaks from macs2 peaks.
        This time we work from merged samples
    Test:
        out/bedtools/intersect_exclude_from_macs2_peaks_in_input_hg38/bedtools/intersect_thymus_stage_exclusive_peaks_hg38/CD34_134.bed
        "Th134_CD34", "Th89_EC", "Th101_LC", "Th125_SP4", "Th125_SP8"
    """
    input:
        bedtools = "opt/bedtools2/bin/bedtools",
        bed_cd34 = "out/macs2/callpeak_format-BAM_gsize-hs_nomodel/ln/alias/experiments/{assembly}_H3K27ac_thymus/CD34_over_input_peaks.narrowPeak",
        bed_ec   = "out/macs2/callpeak_format-BAM_gsize-hs_nomodel/ln/alias/experiments/{assembly}_H3K27ac_thymus/EC_over_input_peaks.narrowPeak",
        bed_lc   = "out/macs2/callpeak_format-BAM_gsize-hs_nomodel/ln/alias/experiments/{assembly}_H3K27ac_thymus/LC_over_input_peaks.narrowPeak",
        bed_sp4  = "out/macs2/callpeak_format-BAM_gsize-hs_nomodel/ln/alias/experiments/{assembly}_H3K27ac_thymus/SP4_over_input_peaks.narrowPeak",
        bed_sp8  = "out/macs2/callpeak_format-BAM_gsize-hs_nomodel/ln/alias/experiments/{assembly}_H3K27ac_thymus/SP8_over_input_peaks.narrowPeak"
    output:
        bed_cd34 = "out/bedtools/intersect_thymus_stage_exclusive_peaks_{assembly}_merged/CD34.bed",
        bed_ec   = "out/bedtools/intersect_thymus_stage_exclusive_peaks_{assembly}_merged/EC.bed",
        bed_lc   = "out/bedtools/intersect_thymus_stage_exclusive_peaks_{assembly}_merged/LC.bed",
        bed_sp4  = "out/bedtools/intersect_thymus_stage_exclusive_peaks_{assembly}_merged/SP4.bed",
        bed_sp8  = "out/bedtools/intersect_thymus_stage_exclusive_peaks_{assembly}_merged/SP8.bed",
    shell:
        """
        {input.bedtools} intersect -a {input.bed_cd34} -b {input.bed_ec} {input.bed_lc} {input.bed_sp4} {input.bed_sp8} -v > {output.bed_cd34}
        {input.bedtools} intersect -a {input.bed_ec} -b {input.bed_cd34} {input.bed_lc} {input.bed_sp4} {input.bed_sp8} -v > {output.bed_ec}
        {input.bedtools} intersect -a {input.bed_lc} -b {input.bed_cd34} {input.bed_ec} {input.bed_sp4} {input.bed_sp8} -v > {output.bed_lc}
        {input.bedtools} intersect -a {input.bed_sp4} -b {input.bed_cd34} {input.bed_ec} {input.bed_lc} {input.bed_sp8} -v > {output.bed_sp4}
        {input.bedtools} intersect -a {input.bed_sp8} -b {input.bed_cd34} {input.bed_ec} {input.bed_lc} {input.bed_sp4} -v > {output.bed_sp8}
        """


rule bedtools_intersect_thymus_stage_common_peaks_hg38_v3:
    """
    Created:
        2018-02-25 15:35:54
    Aim:
        Creating thymus stage common peaks from macs2 peaks.
        V3 because I did not take the samples Salva wanted in the previous intersection.
    Test:
        out/bedtools/intersect_thymus_stage_common_peaks_hg38_v3/multiinter.bed
    """
    input:
        bedtools="opt/bedtools2/bin/bedtools",
        bed_all  ="out/bedtools/intersect_exclude_from_macs2_peaks_in_input_hg38/bedtools/merge/sort/coordinates_bed/cat/hg38-macs2-peaks-H3K27ac-thymus.bed",
        bed_all_subset_3 = "out/bedtools/intersect_exclude_from_macs2_peaks_in_input_hg38/bedtools/merge/sort/coordinates_bed/cat/hg38-macs2-peaks-H3K27ac-thymus-subset-exclusive-3.bed",
        bed_cd34_134="out/sort/coordinates_bed/macs2/callpeak_broad_format-BAM_gsize-hs/ln/alias/experiments/{assembly}_H3K27ac_thymus/Th134_CD34_over_input_peaks.broadPeak",
        bed_ec_89   ="out/sort/coordinates_bed/macs2/callpeak_broad_format-BAM_gsize-hs/ln/alias/experiments/{assembly}_H3K27ac_thymus/Th89_EC_over_input_peaks.broadPeak",
        bed_ec_125  ="out/sort/coordinates_bed/macs2/callpeak_broad_format-BAM_gsize-hs/ln/alias/experiments/{assembly}_H3K27ac_thymus/Th125_EC_over_input_peaks.broadPeak",
        bed_lc_91   ="out/sort/coordinates_bed/macs2/callpeak_broad_format-BAM_gsize-hs/ln/alias/experiments/{assembly}_H3K27ac_thymus/Th91_LC_over_input_peaks.broadPeak",
        bed_lc_101  ="out/sort/coordinates_bed/macs2/callpeak_broad_format-BAM_gsize-hs/ln/alias/experiments/{assembly}_H3K27ac_thymus/Th101_LC_over_input_peaks.broadPeak",
        bed_lc_118  ="out/sort/coordinates_bed/macs2/callpeak_broad_format-BAM_gsize-hs/ln/alias/experiments/{assembly}_H3K27ac_thymus/Th118_LC_over_input_peaks.broadPeak",
        bed_lc_125  ="out/sort/coordinates_bed/macs2/callpeak_broad_format-BAM_gsize-hs/ln/alias/experiments/{assembly}_H3K27ac_thymus/Th125_LC_over_input_peaks.broadPeak",
        bed_sp4_91  ="out/sort/coordinates_bed/macs2/callpeak_broad_format-BAM_gsize-hs/ln/alias/experiments/{assembly}_H3K27ac_thymus/Th91_SP4_over_input_peaks.broadPeak",
        bed_sp4_125 ="out/sort/coordinates_bed/macs2/callpeak_broad_format-BAM_gsize-hs/ln/alias/experiments/{assembly}_H3K27ac_thymus/Th125_SP4_over_input_peaks.broadPeak",
        bed_sp8_125 ="out/sort/coordinates_bed/macs2/callpeak_broad_format-BAM_gsize-hs/ln/alias/experiments/{assembly}_H3K27ac_thymus/Th125_SP8_over_input_peaks.broadPeak"
        ### ADD Peaks from each stage. Need to select which peaks we get from each stage.
    output:
        bed_multiinter="out/bedtools/intersect_thymus_stage_common_peaks_{assembly}_v3/multiinter.bed",
        bed_common="out/bedtools/intersect_thymus_stage_common_peaks_{assembly}_v3/common.bed",
        bed_not_common="out/bedtools/intersect_thymus_stage_common_peaks_{assembly}_v3/not_common.bed",
        bed_not_common_subset_3="out/bedtools/intersect_thymus_stage_common_peaks_{assembly}_v3/not_common_subset_3.bed",
        bed_cd34_134="out/bedtools/intersect_thymus_stage_common_peaks_{assembly}_v3/CD34_134.bed",
        bed_ec_89   ="out/bedtools/intersect_thymus_stage_common_peaks_{assembly}_v3/EC_89.bed",
        ##bed_ec_125  ="out/bedtools/intersect_thymus_stage_exclusive_peaks_{assembly}_v3/EC_125.bed",
        bed_lc_91   ="out/bedtools/intersect_thymus_stage_common_peaks_{assembly}_v3/LC_91.bed",
        ##bed_lc_101  ="out/bedtools/intersect_thymus_stage_exclusive_peaks_{assembly}_v3/LC_101.bed",
        ##bed_lc_118  ="out/bedtools/intersect_thymus_stage_exclusive_peaks_{assembly}_v3/LC_118.bed",
        ##bed_lc_125  ="out/bedtools/intersect_thymus_stage_exclusive_peaks_{assembly}_v3/LC_125.bed",
        bed_sp4_91  ="out/bedtools/intersect_thymus_stage_common_peaks_{assembly}_v3/SP4_91.bed",
        ##bed_sp4_125 ="out/bedtools/intersect_thymus_stage_exclusive_peaks_{assembly}_v3/SP4_125.bed",
        bed_sp8_125 ="out/bedtools/intersect_thymus_stage_common_peaks_{assembly}_v3/SP8_125.bed",
    shell:
        """
        {input.bedtools} multiinter -i {input.bed_cd34_134} {input.bed_ec_89} {input.bed_lc_91} {input.bed_sp4_91} {input.bed_sp8_125} > {output.bed_multiinter}
        awk 'BEGIN {{OFS="\\t"}}; ($4 == '5') {{print $0}}' {output.bed_multiinter} > {output.bed_common}
        {input.bedtools} intersect -a {input.bed_all} -b {output.bed_common} -v > {output.bed_not_common}
        {input.bedtools} intersect -a {input.bed_all_subset_3} -b {output.bed_common} -v > {output.bed_not_common_subset_3}

        {input.bedtools} intersect -a {input.bed_cd34_134} -b {output.bed_common} -v > {output.bed_cd34_134}
        {input.bedtools} intersect -a {input.bed_ec_89} -b {output.bed_common} -v > {output.bed_ec_89}
        {input.bedtools} intersect -a {input.bed_lc_91} -b {output.bed_common} -v > {output.bed_lc_91}
        {input.bedtools} intersect -a {input.bed_sp4_91} -b {output.bed_common} -v > {output.bed_sp4_91}
        {input.bedtools} intersect -a {input.bed_sp8_125} -b {output.bed_common} -v > {output.bed_sp8_125}
        """

rule bedtools_intersect_thymus_stage_common_peaks_hg38_merged:
    """
    Created:
        2018-03-16 19:44:21
    Aim:
        Creating thymus stage common peaks from macs2 peaks.
        Here working with merged samples.
    Test:
        out/bedtools/intersect_thymus_stage_common_peaks_hg38_v3/multiinter.bed
    """
    input:
        bedtools="opt/bedtools2/bin/bedtools",
        bed_all  ="out/bedtools/intersect_exclude_from_macs2_peaks_in_input_hg38/bedtools/merge/sort/coordinates_bed/cat/hg38-macs2-peaks-H3K27ac-thymus-merged.bed",
        bed_cd34 = "out/sort/coordinates_bed/macs2/callpeak_format-BAM_gsize-hs_nomodel/ln/alias/experiments/{assembly}_H3K27ac_thymus/CD34_over_input_peaks.narrowPeak",
        bed_ec   = "out/sort/coordinates_bed/macs2/callpeak_format-BAM_gsize-hs_nomodel/ln/alias/experiments/{assembly}_H3K27ac_thymus/EC_over_input_peaks.narrowPeak",
        bed_lc   = "out/sort/coordinates_bed/macs2/callpeak_format-BAM_gsize-hs_nomodel/ln/alias/experiments/{assembly}_H3K27ac_thymus/LC_over_input_peaks.narrowPeak",
        bed_sp4  = "out/sort/coordinates_bed/macs2/callpeak_format-BAM_gsize-hs_nomodel/ln/alias/experiments/{assembly}_H3K27ac_thymus/SP4_over_input_peaks.narrowPeak",
        bed_sp8  = "out/sort/coordinates_bed/macs2/callpeak_format-BAM_gsize-hs_nomodel/ln/alias/experiments/{assembly}_H3K27ac_thymus/SP8_over_input_peaks.narrowPeak"
    output:
        bed_multiinter = "out/bedtools/intersect_thymus_stage_common_peaks_{assembly}_merged/multiinter.bed",
        bed_common     = "out/bedtools/intersect_thymus_stage_common_peaks_{assembly}_merged/common.bed",
        bed_not_common = "out/bedtools/intersect_thymus_stage_common_peaks_{assembly}_merged/not_common.bed",
        bed_cd34       = "out/bedtools/intersect_thymus_stage_common_peaks_{assembly}_merged/CD34.bed",
        bed_ec         = "out/bedtools/intersect_thymus_stage_common_peaks_{assembly}_merged/EC.bed",
        bed_lc         = "out/bedtools/intersect_thymus_stage_common_peaks_{assembly}_merged/LC.bed",
        bed_sp4        = "out/bedtools/intersect_thymus_stage_common_peaks_{assembly}_merged/SP4.bed",
        bed_sp8        = "out/bedtools/intersect_thymus_stage_common_peaks_{assembly}_merged/SP8.bed",
    shell:
        """
        {input.bedtools} multiinter -i {input.bed_cd34} {input.bed_ec} {input.bed_lc} {input.bed_sp4} {input.bed_sp8} > {output.bed_multiinter}
        awk 'BEGIN {{OFS="\\t"}}; ($4 == '5') {{print $0}}' {output.bed_multiinter} > {output.bed_common}
        {input.bedtools} intersect -a {input.bed_all} -b {output.bed_common} -v > {output.bed_not_common}

        {input.bedtools} intersect -a {input.bed_cd34} -b {output.bed_common} -v > {output.bed_cd34}
        {input.bedtools} intersect -a {input.bed_ec} -b {output.bed_common} -v > {output.bed_ec}
        {input.bedtools} intersect -a {input.bed_lc} -b {output.bed_common} -v > {output.bed_lc}
        {input.bedtools} intersect -a {input.bed_sp4} -b {output.bed_common} -v > {output.bed_sp4}
        {input.bedtools} intersect -a {input.bed_sp8} -b {output.bed_common} -v > {output.bed_sp8}
        """


rule test_dynamic:
    input:
        bed_combination= dynamic("out/bedtools/multiinter_thymus_peaks_hg38_merged/{combination}.bed")

rule bedtools_multiinter_thymus_peaks_hg38_merged_dynamic:
    """
    Created:
        2018-03-19 10:46:58
    Aim:
        Creating thymus stage common peaks from macs2 peaks.
        Here working with merged samples.
        Dynamic working but not completely satisfying.
    Test:
        out/bedtools/multiinter_thymus_peaks_hg38_merged/multiinter.bed
    """
    input:
        bedtools="opt/bedtools2/bin/bedtools",
        bed_all  ="out/bedtools/intersect_exclude_from_macs2_peaks_in_input_hg38/bedtools/merge/sort/coordinates_bed/cat/hg38-macs2-peaks-H3K27ac-thymus-merged.bed",
        bed_cd34 = "out/sort/coordinates_bed/macs2/callpeak_format-BAM_gsize-hs_nomodel/ln/alias/experiments/hg38_H3K27ac_thymus/CD34_over_input_peaks.narrowPeak",
        bed_ec   = "out/sort/coordinates_bed/macs2/callpeak_format-BAM_gsize-hs_nomodel/ln/alias/experiments/hg38_H3K27ac_thymus/EC_over_input_peaks.narrowPeak",
        bed_lc   = "out/sort/coordinates_bed/macs2/callpeak_format-BAM_gsize-hs_nomodel/ln/alias/experiments/hg38_H3K27ac_thymus/LC_over_input_peaks.narrowPeak",
        bed_sp4  = "out/sort/coordinates_bed/macs2/callpeak_format-BAM_gsize-hs_nomodel/ln/alias/experiments/hg38_H3K27ac_thymus/SP4_over_input_peaks.narrowPeak",
        bed_sp8  = "out/sort/coordinates_bed/macs2/callpeak_format-BAM_gsize-hs_nomodel/ln/alias/experiments/hg38_H3K27ac_thymus/SP8_over_input_peaks.narrowPeak"
    output:
        bed_combination= dynamic("out/bedtools/multiinter_thymus_peaks_hg38_merged/{combination}.bed")
    params:
        outdir = "out/bedtools/multiinter_thymus_peaks_hg38_merged",
        bed_multiinter = "out/bedtools/multiinter_thymus_peaks_hg38_merged/multiinter.bed",
        bed_all        = "out/bedtools/multiinter_thymus_peaks_hg38_merged/all_f0.5.bed",
        #bed_not_common = "out/bedtools/multiinter_thymus_peaks_hg38_merged/not_common.bed",
        #bed_cd34       = "out/bedtools/multiinter_thymus_peaks_hg38_merged/CD34.bed",
        #bed_ec         = "out/bedtools/multiinter_thymus_peaks_hg38_merged/EC.bed",
        #bed_lc         = "out/bedtools/multiinter_thymus_peaks_hg38_merged/LC.bed",
        #bed_sp4        = "out/bedtools/multiinter_thymus_peaks_hg38_merged/SP4.bed",
        #bed_sp8        = "out/bedtools/multiinter_thymus_peaks_hg38_merged/SP8.bed",
    shell:
        """
        {input.bedtools} multiinter -i {input.bed_cd34} {input.bed_ec} {input.bed_lc} {input.bed_sp4} {input.bed_sp8} > {params.bed_multiinter}
        {input.bedtools} intersect -a {input.bed_all} -b {params.bed_multiinter} -wa -wb -f 0.5 > {params.bed_all}

        for CLASS in `cut -f8 {params.bed_all} | sort -u`
        do
            echo $CLASS
            awk -v CLASS="$CLASS" 'BEGIN {{OFS="\\t"}}; ($8 == CLASS) {{print $0}}' {params.bed_all} > {params.outdir}/$CLASS.bed
        done

        rm -f {params.bed_multiinter} {params.bed_all}
        """



rule bedtools_intersect_bs_hypometh_thymus_merge_and_multiinter:
    input:
        a = "out/bedtools/merge/bedtools/multiinter_bed-hg38-bs-hypometh-thymus/multiinter.bed",
        b = "out/bedtools/multiinter_bed-hg38-bs-hypometh-thymus/multiinter.bed"
    output:
        bed = "out/bedtools/intersect_bs_hypometh_thymus_merge_and_multiinter/test.bed"
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        bedtools intersect -wa -wb -a {input.a} -b {input.b} > {output.bed}
        """
