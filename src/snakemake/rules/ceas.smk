rule ceas:
    """
    Created:
        2017-02-21 15:33:10
    Aim:
        Standard CEAS analysis for one sample
    Note:
        Unfinished
    Test:
        out/ceas/refgene-mm9_bed-mm9-test-macs2-peaks/danpos/dtriple/samtools/index/samtools/sort/samtools/view_sam_to_bam_-q_30/bowtie2/se_mm9/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR3126243.R
        Error: I produce SRR3126243_over_SRR3126242_peaks.R

    http://liulab.dfci.harvard.edu/CEAS/src/mm9.refGene.gz
        out/ceas/refgene-mm9_bed-mm9-test-macs2-peaks/danpos/dtriple/samtools/index/samtools/sort/samtools/view_sam_to_bam_-q_30/bowtie2/se_mm9/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR3126243.R
    """
    input:
        average_profiling        = "out/{filler}.wig",
        chip_region_annotation   = lambda wildcards: eval(mwconf['ids'][wildcards.bed_id]),
        gene_centered_annotation = lambda wildcards: eval(mwconf['ids'][wildcards.refgene_id])
    output:
        rcode = "out/ceas/{refgene_id}_{bed_id}/{filler}.R",
        #pdf   = "out/ceas/feature/{index}/{rep_pos}/{exp}/{sample}/{selection}/{selection}.pdf"
    conda:
        "../envs/ceas.yaml"
    threads:
        1
    shell:
        """
        WDIR=`pwd`
        OUTDIR=`dirname {output.rcode}`
        rm -rf $OUTDIR
        mkdir -p $OUTDIR
        cd $OUTDIR
        
        ceas \
            -g $WDIR/{input.gene_centered_annotation} \
            -b $WDIR/{input.chip_region_annotation} \
            -w $WDIR/{input.average_profiling}
        
        # DEBUG
        #cp -r $WDIR/outdir $WDIR/tmp/debug/
        """

rule for_article_4a_diff_clust_mm10_PSKs_SC_WT_fuzz_lt_40_100bp_kmeans20:
    """
    Created: 2016-05-30 16h15
    
    Visually, we identify clusters 6,8,10,12,16,18 where signal is higher in WT. In the same way, we identify clusters 5,9,11 where signal is higher in KO. Other clusters have signal similar in both samples.
    """
    input:
        bed="result/deepTools/plotHeatmap/referencePoint/center/bamCoverage/mm10/for_article_4a_features/top_pos/merge_run113_run125_run126_PSK-SC-WT_fuzz_lt_40/100bp/kmeans20/yMax200/all.bed"
    output:
        wt_lt_ko="annotation/processed/feature/mm10/diff_clust/PSK-SC-WT_fuzz_lt_40_100bp_kmeans20/WT_lt_KO.bed",
        wt_gt_ko="annotation/processed/feature/mm10/diff_clust/PSK-SC-WT_fuzz_lt_40_100bp_kmeans20/WT_gt_KO.bed",
        wt_eq_ko="annotation/processed/feature/mm10/diff_clust/PSK-SC-WT_fuzz_lt_40_100bp_kmeans20/WT_eq_KO.bed"
    shell:"""
    awk 'BEGIN{{FS=OFS="\t"}}$7 ~ "cluster [6 8 10 12 16 18]"{{print $1,$2,$3}}' {input.bed} > {output.wt_gt_ko}
    awk 'BEGIN{{FS=OFS="\t"}}$7 ~ "cluster [5 9 11]"{{print $1,$2,$3}}' {input.bed} > {output.wt_lt_ko}
    awk 'BEGIN{{FS=OFS="\t"}}$7 ~ "cluster [1 2 3 4 5 7 13 14 15 17 19 20]"{{print $1,$2,$3}}' {input.bed} > {output.wt_eq_ko}
    """

rule for_article_4a_diff_clust_mm9_PSKs_SC_WT_fuzz_lt_40_100bp_kmeans20:
    """
    Created: 2016-05-30 16h15
    
    Visually, we identify clusters 6,9,12,15,17 where signal is higher in WT. In the same way, we identify clusters 7,8,10 where signal is higher in KO. Other clusters have signal similar in both samples.
    """
    input:
        bed="result/deepTools/plotHeatmap/referencePoint/center/bamCoverage/mm9/for_article_4a_features/top_pos/merge_run113_run125_run126_PSK-SC-WT_fuzz_lt_40/100bp/kmeans20/yMax200/all.bed"
    output:
        wt_lt_ko="annotation/processed/feature/mm9/diff_clust/PSK-SC-WT_fuzz_lt_40_100bp_kmeans20/WT_lt_KO.bed",
        wt_gt_ko="annotation/processed/feature/mm9/diff_clust/PSK-SC-WT_fuzz_lt_40_100bp_kmeans20/WT_gt_KO.bed",
        wt_eq_ko="annotation/processed/feature/mm9/diff_clust/PSK-SC-WT_fuzz_lt_40_100bp_kmeans20/WT_eq_KO.bed"
    shell:"""
    awk 'BEGIN{{FS=OFS="\t"}}$7 ~ "cluster [6 9 12 15 17]"{{print $1,$2,$3}}' {input.bed} > {output.wt_gt_ko}
    awk 'BEGIN{{FS=OFS="\t"}}$7 ~ "cluster [7 8 10]"{{print $1,$2,$3}}' {input.bed} > {output.wt_lt_ko}
    awk 'BEGIN{{FS=OFS="\t"}}$7 ~ "cluster [1 2 3 4 5 11 13 14 16 18 19 20]"{{print $1,$2,$3}}' {input.bed} > {output.wt_eq_ko}
    """


rule for_article_4a_diff_clust_mm9_PSKs_SC_WT_fuzz_lt_45_100bp_kmeans20:
    """
    Created: 2016-05-30 16h15
    Modified: 2016-05-30 17h46 - Actually it is better to work with mm9 so I can use the builin refgene for CEAS.
    
    Working with 45 as fuzziness threshold and 20 k-mean classes on mm9 assembly, we get figure~\ref{fig:profile-and-heatmap-for-psk-sc-wt-and-psk-sc-ko-around-nucleosomes-with-fuzziness-lower-than-45-kmeans-20-mm9}. Visually, we identify clusters 8,10,13,14,18 where signal is higher in WT. In the same way, we identify clusters 9,11,12 where signal is higher in KO. Other clusters have signal similar in both samples.
    
    """
    input:
        bed="result/deepTools/plotHeatmap/referencePoint/center/bamCoverage/mm9/for_article_4a_features/top_pos/merge_run113_run125_run126_PSK-SC-WT_fuzz_lt_45/100bp/kmeans20/yMax200/all.bed"
    output:
        wt_lt_ko="annotation/processed/feature/mm9/diff_clust/PSK-SC-WT_fuzz_lt_45_100bp_kmeans20/WT_lt_KO.bed",
        wt_gt_ko="annotation/processed/feature/mm9/diff_clust/PSK-SC-WT_fuzz_lt_45_100bp_kmeans20/WT_gt_KO.bed",
        wt_eq_ko="annotation/processed/feature/mm9/diff_clust/PSK-SC-WT_fuzz_lt_45_100bp_kmeans20/WT_eq_KO.bed"
    shell:"""
    awk 'BEGIN{{FS=OFS="\t"}}$7 ~ "cluster [8 10 13 14 18]"{{print $1,$2,$3}}' {input.bed} > {output.wt_gt_ko}
    awk 'BEGIN{{FS=OFS="\t"}}$7 ~ "cluster [9 11 12]"{{print $1,$2,$3}}' {input.bed} > {output.wt_lt_ko}
    awk 'BEGIN{{FS=OFS="\t"}}$7 ~ "cluster [1 2 3 4 5 6 7 10 14 15 16 17 19 20]"{{print $1,$2,$3}}' {input.bed} > {output.wt_eq_ko}
    """

rule for_article_4d_diff_clust_mm9_PSKs_SC_WT_fuzz_lt_45_100bp_kmeans20:
    """
    Created: 2016-05-30 16h15
    Modified: 2016-05-30 17h46 - Actually it is better to work with mm9 so I can use the builin refgene for CEAS.
    
    Working with 45 as fuzziness threshold and 20 k-mean classes on mm9 assembly, we get figure~\ref{fig:profile-and-heatmap-for-psk-sc-wt-psk-sc-ko-and-its-two-subpopulations-around-nucleosomes-with-fuzziness-lower-than-45-kmeans-20-mm9}. Visually, we identify clusters 5,6,12,13,14,18 where signal is higher in 47-59bp than in 59-71bp. In the same way, we identify clusters 7,9,11,15,16,17 where signal is higher in 59-71bp. Other clusters have signal similar in both samples.
    
    """
    input:
        bed="result/deepTools/plotHeatmap/referencePoint/center/bamCoverage/mm9/for_article_4d_features/top_pos/merge_run113_run125_run126_PSK-SC-WT_fuzz_lt_45/100bp/kmeans20/yMax200/all.bed"
    output:
        wt_lt_ko="annotation/processed/feature/mm9/diff_clust/PSK-SC-WT_fuzz_lt_45_100bp_kmeans20_4d/47-59_lt_59-71.bed",
        wt_gt_ko="annotation/processed/feature/mm9/diff_clust/PSK-SC-WT_fuzz_lt_45_100bp_kmeans20_4d/47-59_gt_59-71.bed",
        wt_eq_ko="annotation/processed/feature/mm9/diff_clust/PSK-SC-WT_fuzz_lt_45_100bp_kmeans20_4d/47-59_eq_59-71.bed"
    shell:"""
    awk 'BEGIN{{FS=OFS="\t"}}$7 ~ "cluster [5 6 12 13 14 18]"{{print $1,$2,$3}}' {input.bed} > {output.wt_gt_ko}
    awk 'BEGIN{{FS=OFS="\t"}}$7 ~ "cluster [7 9 11 15 16 17]"{{print $1,$2,$3}}' {input.bed} > {output.wt_lt_ko}
    awk 'BEGIN{{FS=OFS="\t"}}$7 ~ "cluster [1 2 3 4 8 10 19 20]"{{print $1,$2,$3}}' {input.bed} > {output.wt_eq_ko}
    """

rule ceas_on_diff_clust:
    """
    Created: 2016-05-23 18h30 from old ceas_danpos_new rule
    
    The aim is to use CEAS (Cis-regulatory Element Annotation System) to assess enrichment in some specific called nucleosomes instead of the whole sample.
    
    "result/ceas/feature/{index}/{rep_pos}/{exp}/{sample}/{selection}/{sample}_{selection}.pdf"
    
    {rep_pos}="rep_pos" right now but could be extended to other features.
    """
    input:
        ceas="opt/anaconda2/bin/ceas",
        chip_region_annotation="annotation/processed/feature/{index}/diff_clust/{diff_clust_sett}/{selection}.bed",
        gene_centered_annotation="annotation/input/{index}.refGene"
    output:
        outdir="result/ceas/feature/{index}/diff_clust/{diff_clust_sett}/{selection}/",
        rcode="result/ceas/feature/{index}/diff_clust/{diff_clust_sett}/{selection}/{selection}.R",
        pdf="result/ceas/feature/{index}/diff_clust/{diff_clust_sett}/{selection}/{selection}.pdf"
    threads: 1
    shell:"""
    WDIR=`pwd`
    rm -rf {output.outdir}
    mkdir -p {output.outdir}
    cd {output.outdir}
    
    $WDIR/{input.ceas} -g $WDIR/{input.gene_centered_annotation} \
    -b $WDIR/{input.chip_region_annotation} 
    """

rule ceas_on_diff_class_fig3_alt:
    """
    Created: 2016-06-23 11h39
    Using CEAS to see if the small structures present in KO and not WT (and vice-versa) have an enrichment bias in some type of regions. 
    """
    input:
        ceas="opt/anaconda2/bin/ceas",
        chip_region_annotation="result/for_article/questions_v5/fig3_alt/pos_PSK-SC-WT_{lt_gt}_KO.bed",
        gene_centered_annotation="annotation/input/mm9.refGene"
    output:
        outdir="result/ceas/feature/mm9/diff_class_fig3_alt/pos_PSK-SC-WT_{lt_gt}_KO/",
        rcode="result/ceas/feature/mm9/diff_class_fig3_alt/pos_PSK-SC-WT_{lt_gt}_KO/pos_PSK-SC-WT_{lt_gt}_KO.R",
        pdf="result/ceas/feature/mm9/diff_class_fig3_alt/pos_PSK-SC-WT_{lt_gt}_KO/pos_PSK-SC-WT_{lt_gt}_KO.pdf"
    threads: 1
    shell:"""
    WDIR=`pwd`
    rm -rf {output.outdir}
    mkdir -p {output.outdir}
    cd {output.outdir}
    
    $WDIR/{input.ceas} -g $WDIR/{input.gene_centered_annotation} \
    -b $WDIR/{input.chip_region_annotation} 
    """


rule ceas_macs:
    """
    Cis-regulatory Element Annotation System
    NOTE: actually part of CEAS analysis with wig file does not work... Check that later or just use other tools like gtftoolkit tss_plot to do the same.
    """
    input:
        chip_region_annotation="out/macs14/{id}/{id}_peaks.bed",
        average_profiling="out/bw/bamCoverage/{id}.bw",
        gene_centered_annotation="annotation/mm9.refGene"
    output: "result/ceas/{id}/done"
    shell:"""
    mkdir result/ceas result/ceas/{wildcards.id}
    cd result/ceas/{wildcards.id}
    # /cobelix/Charbonnier/.local/bin/ceas
    ceas -g ../../../{input.gene_centered_annotation} -b ../../../{input.chip_region_annotation} -w ../../../{input.average_profiling}
    
    Rscript {wildcards.id}_peaks.R
    
    touch ../../../{output}
    """

# 2016-07-22 11h54 - Moved this part of the code in install_opt.rules. Remove it here if no problem seen later. 
#rule install_ceas:
#	"""
#	Created: 2016-04-12 15h18
#	Not used right now.
#	I wanted it to use build_genomeBG for mm10 but it needs wig not available for mm10...
#	
#	Maybe obsolete as it is easier to install with bioconda.
#	Note, this install seems to put lib in /cobelix/Charbonnier/.local/lib/python2.7/site-packages/CEAS/inout.py:64: UserWarning: sqlite3 is used instead of MySQLdb because MySQLdb is not installed
#	"""
#	input:
#	output: ceas="opt/anaconda2/bin/ceas"
#	shell:"""
#	cd opt
#	wget http://liulab.dfci.harvard.edu/CEAS/src/CEAS-Package-1.0.2.tar.gz
#	tar xvf CEAS-Package-1.0.2.tar.gz 
#	cd CEAS-Package-1.0.2
#	../anaconda2/bin/python setup.py install
#	
#	"""
#
#rule install_ceas_alt_by_conda:
#	"""
#	Created: 2016-05-25 11h02
#	Created to check if a standard install by conda put the library in the project directory.
#	"""


rule ceas_on_features:
    """
    Created: 2016-05-23 18h30 from old ceas_danpos_new rule
    
    The aim is to use CEAS (Cis-regulatory Element Annotation System) to assess enrichment in some specific called nucleosomes instead of the whole sample.
    
    "result/ceas/feature/{index}/{rep_pos}/{exp}/{sample}/{selection}/{sample}_{selection}.pdf"
    
    {rep_pos}="rep_pos" right now but could be extended to other features.
    """
    input:
        ceas="opt/anaconda2/bin/ceas",
        chip_region_annotation="annotation/processed/feature/{index}/{rep_pos}/{exp}/{sample}/{selection}.bed",
        average_profiling="out/danpos/dtriple/{index}/{exp}/{sample}/pooled/data_processed_bam_{index}_{exp}_{sample}.smooth.wig",
        gene_centered_annotation="annotation/input/{index}.refGene"
    output:
        outdir="result/ceas/feature/{index}/{rep_pos}/{exp}/{sample}/{selection}/",
        rcode="result/ceas/feature/{index}/{rep_pos}/{exp}/{sample}/{selection}/{selection}.R",
        pdf="result/ceas/feature/{index}/{rep_pos}/{exp}/{sample}/{selection}/{selection}.pdf"
    threads: 1
    shell:"""
    WDIR=`pwd`
    rm -rf {output.outdir}
    mkdir -p {output.outdir}
    cd {output.outdir}
    
    $WDIR/{input.ceas} -g $WDIR/{input.gene_centered_annotation} \
    -b $WDIR/{input.chip_region_annotation} \
    -w $WDIR/{input.average_profiling}
    
    # DEBUG
    cp -r $WDIR/{output.outdir} $WDIR/tmp/debug/
    """

rule ceas_danpos_new:
    """
    Created: 2016-04-12 15h10 from old ceas_danpos rule
    
    Cis-regulatory Element Annotation System
    NOTE: actually part of CEAS analysis with wig file does not work... Check that later or just use other tools like gtftoolkit tss_plot to do the same.
    # TRY TO ADD THAT LATER TEST pdf="result/ceas/danpos/{sample}/{ppr}/{sample}_{ppr}.pdf"
    expand("result/ceas/danpos/dtriple/{index}/{exp}/{sample}/{ppr}/{sample}_{ppr}.R", index="mm9", exp="merge_run113_run119_run124", sample="MNS-R-WT", ppr="positions")
    """
    input:
        ceas="opt/miniconda/envs/py27/bin/ceas",
        chip_region_annotation="out/danpos/dtriple/{index}/{exp}/{sample}/pooled/data_processed_bam_{index}_{exp}_{sample}.smooth.{ppr}.xls",
        average_profiling="out/danpos/dtriple/{index}/{exp}/{sample}/pooled/data_processed_bam_{index}_{exp}_{sample}.smooth.wig",
        gene_centered_annotation="annotation/input/{index}.refGene"
    output:
        outdir="result/ceas/danpos/dtriple/{index}/{exp}/{sample}/{ppr}/",
        rcode="result/ceas/danpos/dtriple/{index}/{exp}/{sample}/{ppr}/{sample}_{ppr}.R"
    shell:"""
    WDIR=`pwd`
    rm -rf {output.outdir}
    mkdir -p {output.outdir}
    cd {output.outdir}
    
    tail -n +2 $WDIR/{input.chip_region_annotation} | cut -f1-3 > {wildcards.sample}_{wildcards.ppr}.bed
    
    $WDIR/{input.ceas} -g $WDIR/{input.gene_centered_annotation} \
    -b {wildcards.sample}_{wildcards.ppr}.bed \
    -w $WDIR/{input.average_profiling}
    
    #Rscript *.R
    
    rm -f {wildcards.sample}_{wildcards.ppr}.bed
    """

rule bedtools_define_distal_intergenic_for_ceas:
    """
    Created: 2017-01-08 13h06
    """
    input:
        gene_bed="annotation/processed/feature/mm9/ensGene.bed",
        mappable_genome_bed="out/awk/get_non_unmapable/ucsc/bigWigToBedGraph/rsync/ucsc/gbdb/mm9/bbi/crgMapabilityAlign75mer.bedGraph",
        chrominfo="out/mysql/ucsc/mm9/chromInfo_main.tsv",
        bedtools="opt/bedtools2/bin/bedtools"
    output:
        distal_intergenic_bed="out/bedtools/define_distal_intergenic_for_ceas/distal_intergenic.bed",
        merge_distal_intergenic_bed="out/bedtools/define_distal_intergenic_for_ceas/merge_distal_intergenic.bed",
        slop_gene_bed="out/bedtools/define_distal_intergenic_for_ceas/slop_gene_3000b.bed"
    shell:"""
    {input.bedtools} slop -i {input.gene_bed} -g {input.chrominfo} -b 3000 > {output.slop_gene_bed}

    {input.bedtools} subtract -a {input.mappable_genome_bed} -b {output.slop_gene_bed} > {output.distal_intergenic_bed}
    
    {input.bedtools} merge -i {output.distal_intergenic_bed} > {output.merge_distal_intergenic_bed}
    """


rule tmp_ceas_test_ebed:
    input:
        expand("result/ceas/test_ebed/danpos/dtriple/mm9/merge_run113_run119_run124/{id}/positions/ceas.R",
            id=["MNS-P-WT", "MNS-P-KO", "MNS-R-WT", "MNS-R-KO"]),
        expand("result/ceas/test_ebed/danpos/dtriple/mm9/merge_run113_run125_run126/{id}/positions/ceas.R",
            id=["MNS-SC-WT", "MNS-SC-KO", "PSK-SC-WT", "PSK-SC-KO"])

rule test_extract_features_from_ceas_features:
    input:
        sqlite3="",
        refgene=""
    output:
    shell:"""
    """

rule ceas_test_ebed:
    """
    Created: 2017-01-08 11h24

    Test:
    "result/ceas/test_ebed/danpos/dtriple/mm9/merge_run113_run119_run124/MNS-P-WT/positions/ceas.R"
    """
    input:
        ceas="opt/miniconda/envs/py27/bin/ceas",
        chip_region_annotation="out/danpos/dtriple/{index}/{exp}/{sample}/pooled/data_processed_bam_{index}_{exp}_{sample}.smooth.{ppr}.xls",
        gene_centered_annotation="annotation/input/mm9.refGene",
        ebed="out/bedtools/define_distal_intergenic_for_ceas/merge_distal_intergenic.bed"
    output:
        rcode="result/ceas/test_ebed/danpos/dtriple/{index}/{exp}/{sample}/{ppr}/ceas.R"
    params:
        outdir="result/ceas/test_ebed/danpos/dtriple/{index}/{exp}/{sample}/{ppr}/"
    shell:"""
    WDIR=`pwd`
    rm -rf {params.outdir}
    mkdir -p {params.outdir}
    cd {params.outdir}
    
    tail -n +2 $WDIR/{input.chip_region_annotation} | cut  -f1-3 > $WDIR/{params.outdir}/chip_region_annotation.bed
    # With true bed, CEAS complains that it can not convert string '.' to float, so we remove extra columns.
    cut -f1-3 $WDIR/{input.ebed} > $WDIR/{params.outdir}/ebed.bed
       
 
    $WDIR/{input.ceas} \
        --gt $WDIR/{input.gene_centered_annotation} \
        --bed $WDIR/{params.outdir}/chip_region_annotation.bed \
        --ebed $WDIR/{params.outdir}/ebed.bed\
        --name ceas
    """

rule manual_ceas_r_barplot_data_for_article:
    input:
        #ceas_r="result/ceas/{filler}.R",
        Rscript="opt/miniconda/envs/r/bin/Rscript",
        code="code/r/script/manual_ceas_r_barplot_data_for_article.R"
    output:
        #tsv="result/ceas/{filler}_barplot.tsv"
        "test_merge_barplot_promoter.pdf"
    # Conda envs not working because r-base needs to be rebuild with conda-build 2.0
    #conda:
    #    "code/snakemake/envs/r.yaml"
    shell:
        """
        {input.Rscript} {input.code}
        """

rule custom_ceas_barplot_for_article:
    """
    Created: 2016-12-1 9h38 
    """
    input:
        tsv="result/ceas/parse_ceas_barplot_data.tsv",
        tsv_lab="result/ceas/parse_ceas_barplot_data_lab.tsv",
        code="code/r/script/custom_ceas_barplot.R",
        rscript="opt/miniconda/envs/r/bin/Rscript"
    output:
        pdf="result/ceas/custom_ceas_barplot_for_article.pdf"
    shell:"""
    {input.rscript} {input.code} -i {input.tsv} -l {input.tsv_lab} -o {output.pdf}
    """

rule parse_ceas_barplot_data:
    """
    Created: 2016-11-29 9h37 - Parsing CEAS barplot data included in multiple standard *.R files produced by CEAS.
    Should only work for mm9 standard reference file as there are hardcoded line selection
    Produce two CSV, one containing 
    """
    input:
        ceas_r=expand(
            "result/ceas/danpos/dtriple/{index}/{exp}/{sample}/{ppr}/{sample}_{ppr}.R",
            zip,
            index=["mm9"]*8,
            exp=["merge_run113_run119_run124"]*4 + ["merge_run113_run125_run126"]*4,
            sample=["MNS-P-WT", "MNS-P-KO", "MNS-R-WT", "MNS-R-KO", "MNS-SC-WT", "MNS-SC-KO", "PSK-SC-WT", "PSK-SC-KO"],
            ppr=["positions"]*8
            )
    output:
        csv="result/ceas/parse_ceas_barplot_data.csv",
        tsv="result/ceas/parse_ceas_barplot_data.tsv",
        csv_lab="result/ceas/parse_ceas_barplot_data_lab.csv",
        tsv_lab="result/ceas/parse_ceas_barplot_data_lab.tsv"
    params:
        outdir="result/ceas/parse_ceas_barplot_data"
    shell:"""
    mkdir -p {params.outdir}

    # Extracting percentages and pvalues from CEAS R output and merging them for multiple samples.
    # This is done with the assumption that lines in CEAS R output does not change for different samples. Let me know if it is not the case.
    
    ### Genome reference
    #### Only retrieve information on the first sample as they are the same.
    FIRST_SAMPLE=`echo "{input.ceas_r}" | cut -f1 --delimiter=" "`

    # Feature names
    ## Chromosomes
    sed -n '34,34 p' $FIRST_SAMPLE | sed 's/^names=c(/Feature,/g' | sed 's/)$//g' | sed 's/,"/,"chr/g' > {params.outdir}/chr.csv
    ## Promoters
    sed -n '48,48 p' $FIRST_SAMPLE | sed 's/^names=c(/Feature,/g' | sed 's/)$//g' | sed 's/,"/,"prom/g' > {params.outdir}/prom.csv
    ## Bidirectional promoters 
    sed -n '59,59 p' $FIRST_SAMPLE | sed 's/^names=c(/Feature,/g' | sed 's/)$//g' | sed 's/,"/,"bidirprom/g' > {params.outdir}/bidirprom.csv
    ## Downstream 
    sed -n '69,69 p' $FIRST_SAMPLE | sed 's/^names=c(/Feature,/g' | sed 's/)$//g' | sed 's/,"/,"downstream/g' > {params.outdir}/downstream.csv
    ## Genes 
    sed -n '80,80 p' $FIRST_SAMPLE | sed 's/^names=c(/Feature,/g' | sed 's/)$//g' | sed 's/,"/,"gene/g' > {params.outdir}/gene.csv

    # Get same header for label table
    cp {params.outdir}/chr.csv {params.outdir}/chr_lab.csv
    cp {params.outdir}/prom.csv {params.outdir}/prom_lab.csv
    cp {params.outdir}/bidirprom.csv {params.outdir}/bidirprom_lab.csv
    cp {params.outdir}/downstream.csv {params.outdir}/downstream_lab.csv
    cp {params.outdir}/gene.csv {params.outdir}/gene_lab.csv

    # Percentage values for reference genome
    ## Chromosomes
    sed -n '36,36 p' $FIRST_SAMPLE | sed 's/^text(x=c(/Genome,/g' | sed 's/),y=mp.*$//g' | sed 's/ //g' >> {params.outdir}/chr.csv
    ##  Promoters (remember to escape [] in sed)
    sed -n '50,50 p' $FIRST_SAMPLE | sed 's/^text(x=mp\[1,\],y=c(/Genome,/g' | sed 's/),label=c(.*$//g' | sed 's/ //g' >> {params.outdir}/prom.csv
    ## Bidirectional promoters
    sed -n '61,61 p' $FIRST_SAMPLE | sed 's/^text(x=mp\[1,\],y=c(/Genome,/g' | sed 's/),label=c(.*$//g' | sed 's/ //g' >> {params.outdir}/bidirprom.csv
    ## Downstream
    sed -n '71,71 p' $FIRST_SAMPLE | sed 's/^text(x=mp\[1,\],y=c(/Genome,/g' | sed 's/),label=c(.*$//g' | sed 's/ //g' >> {params.outdir}/downstream.csv
    ## Gene
    sed -n '82,82 p' $FIRST_SAMPLE | sed 's/^text(x=mp\[1,\],y=c(/Genome,/g' | sed 's/),label=c(.*$//g' | sed 's/ //g' >> {params.outdir}/gene.csv

    # Labels for reference genome
    sed -n '36,36 p' $FIRST_SAMPLE | sed 's/^.*label=c(/Genome,/g' | sed 's/),pos=4.*$//g' | sed 's/ //g' | sed 's/%/% /g' >> {params.outdir}/chr_lab.csv
    ##  Promoters (remember to escape [] in sed)
    sed -n '50,50 p' $FIRST_SAMPLE | sed 's/^.*label=c(/Genome,/g' | sed 's/),pos=3.*$//g' | sed 's/ //g' | sed 's/%/% /g' >> {params.outdir}/prom_lab.csv
    ## Bidirectional promoters
    sed -n '61,61 p' $FIRST_SAMPLE | sed 's/^.*label=c(/Genome,/g' | sed 's/),pos=3.*$//g' | sed 's/ //g' | sed 's/%/% /g' >> {params.outdir}/bidirprom_lab.csv
    ## Downstream
    sed -n '71,71 p' $FIRST_SAMPLE | sed 's/^.*label=c(/Genome,/g' | sed 's/),pos=3.*$//g' | sed 's/ //g' | sed 's/%/% /g' >> {params.outdir}/downstream_lab.csv
    ## Gene
    sed -n '82,82 p' $FIRST_SAMPLE | sed 's/^.*label=c(/Genome,/g' | sed 's/),pos=3.*$//g' | sed 's/ //g' | sed 's/%/% /g' >> {params.outdir}/gene_lab.csv

    ## Samples:
    for SAMPLE in {input.ceas_r}
    do
        # Get sample name for table
        NAME=`sed -n '18,18 p' $SAMPLE | sed -e 's/pdf("\(.*\)_positions.pdf.*$/\\1/'`
        # Percentage values for sample using NAME
        ## Chromosomes
        sed -n '37,37 p' $SAMPLE | sed "s/^text(x=c(/$NAME,/g" | sed 's/),y=mp.*$//g' | sed 's/ //g' >> {params.outdir}/chr.csv
        ## Promoters (labels are on multiple lines for samples)
        sed -n '51,51 p' $SAMPLE | sed "s/^text(x=mp\[2,\],y=c(/$NAME,/g" | sed 's/),label=.*$//g' | sed 's/ //g' >> {params.outdir}/prom.csv
        ## Bidirectional promoters
        sed -n '62,62 p' $SAMPLE | sed "s/^text(x=mp\[2,\],y=c(/$NAME,/g" | sed 's/),label=.*$//g' | sed 's/ //g' >> {params.outdir}/bidirprom.csv
        ## Downstream
        sed -n '72,72 p' $SAMPLE | sed "s/^text(x=mp\[2,\],y=c(/$NAME,/g" | sed 's/),label=.*$//g' | sed 's/ //g' >> {params.outdir}/downstream.csv
        ## Gene
        sed -n '83,83 p' $SAMPLE | sed "s/^text(x=mp\[2,\],y=c(/$NAME,/g" | sed 's/),label=.*$//g' | sed 's/ //g' >> {params.outdir}/gene.csv

        # Labels for sample using NAME
        ## Chromosomes
        sed -n '37,37 p' $SAMPLE | sed "s/^text(x=c(/$NAME,/g" | sed 's/),y=mp.*$//g' | sed 's/ //g' >> {params.outdir}/chr_lab.csv
        ## Promoters (labels are on multiple lines for samples)
        sed -n '51,54 p' $SAMPLE | tr -d '\n' | sed "s/^.*label=c(/$NAME,/g" | sed -e 's/),pos=.*$/\\n/g' | sed 's/ //g'| sed 's/%/% /g' >> {params.outdir}/prom_lab.csv
        ## Bidirectional promoters
        sed -n '62,64 p' $SAMPLE | tr -d '\n' | sed "s/^.*label=c(/$NAME,/g" | sed -e 's/),pos=.*$/\\n/g' | sed 's/ //g'| sed 's/%/% /g' >> {params.outdir}/bidirprom_lab.csv
        ## Downstream
        sed -n '72,75 p' $SAMPLE | tr -d '\n' | sed "s/^.*label=c(/$NAME,/g" | sed -e 's/),pos=.*$/\\n/g' | sed 's/ //g'| sed 's/%/% /g' >> {params.outdir}/downstream_lab.csv
        ## Gene
        sed -n '83,88 p' $SAMPLE | tr -d '\n' | sed "s/^.*label=c(/$NAME,/g" | sed -e 's/),pos=.*$/\\n/g' | sed 's/ //g'| sed 's/%/% /g' >> {params.outdir}/gene_lab.csv
    done

    paste --delimiters ',' {params.outdir}/chr.csv {params.outdir}/prom.csv {params.outdir}/bidirprom.csv {params.outdir}/downstream.csv {params.outdir}/gene.csv > {output.csv}
    paste --delimiters ',' {params.outdir}/chr_lab.csv {params.outdir}/prom_lab.csv {params.outdir}/bidirprom_lab.csv {params.outdir}/downstream_lab.csv {params.outdir}/gene_lab.csv > {output.csv_lab}

    sed 's/,/\\t/g' {output.csv} > {output.tsv}
    sed 's/,/\\t/g' {output.csv_lab} > {output.tsv_lab}
    """

rule merge_ceas_r_pie_data:
    """
    Created: 2016-07-25 9h55 - Merging CEAS pie data to produce a more compact figure for article.
    """
    input:
        ceas_r=expand("result/ceas/danpos/dtriple/{index}/{exp}/{sample}/{ppr}/{sample}_{ppr}.R", zip, index=["mm9"]*8, exp=["merge_run113_run119_run124"]*4+["merge_run113_run125_run126"]*4, sample=["MNS-P-WT", "MNS-P-KO", "MNS-R-WT", "MNS-R-KO", "MNS-SC-WT", "MNS-SC-KO", "PSK-SC-WT", "PSK-SC-KO"], ppr=["positions"]*8)
    output:
        merge_r="result/ceas/merge_r_pie_data.R",
        csv="result/ceas/merge_r_pie_data.csv",
        tsv="result/ceas/merge_r_pie_data.tsv"
    shell:"""
    # Printing the line containing the data for pie plots.
    ## Genome reference
    ### head -1 to select only one because they are all the same.
    ### sed to change the name of the sample into 'Genome'.
    grep "pie(.*Genome" {input.ceas_r} | head -1 | sed 's/.*_positions.R/Genome_positions.R/' > {output.merge_r}
    ## 
    # Samples. (CEAS call them ChIP by default)
    grep "pie(.*ChIP" {input.ceas_r} >> {output.merge_r}
    
    # Retrieving sample names.
    cat {output.merge_r} | sed 's/\(.*\)_positions.R.*$/\\1/' | awk -F"/" '{{print $(NF)}}' > {output.merge_r}_names.tmp
    
    # Retrieving data from labels in R code.
    cat {output.merge_r} | sed 's/^.*labels=c(\(.*\)),main.*$/\\1/' | sed 's/["% ]//g' > {output.merge_r}_data.tmp
    
    # Printing labels
    echo "Class,Prom <1kb,Prom 1-2kb,Prom 2-3kb,DS <1kb,DS 1-2kb,DS 2-3kb,5'UTR,3'UTR,Cod.exons,Introns,Intergenic" > {output.csv}
    paste --delimiters=',' {output.merge_r}_names.tmp {output.merge_r}_data.tmp >> {output.csv}
    
    sed 's/,/\t/g' {output.csv} > {output.tsv}
    """

rule pie_chart_merge_ceas:
    input:
        code="code/r/pie_chart_merge_ceas.R",
        rscript="opt/miniconda/envs/py35/bin/Rscript",
        tsv="result/ceas/merge_r_pie_data.tsv"
    output:
        pdf="result/ceas/merge_r_pie_data.pdf"
    shell:"""
    {input.rscript} {input.code}
    """

rule merge_ceas_r_output:
    """
    The profiles are the same for positions, peaks and regions as the same wig is used for each.
    
    input commented to work on this rule on pedagogix. Remove comment to restore workflow functionality.
    """
    input:
        indir="result/ceas/danpos"#,\
        #rcode=expand("result/ceas/danpos/{sample}/{ppr_in}/{sample}_{ppr_in}.R", ppr_in=PPR, sample=SAMPLES_MNASE)
    output:
        rcode="result/merge_ceas_r_output/{ppr}_code.r"
    shell:"""
    OUTDIR="result/merge_ceas_r_output"
    SAMPLES="{input.indir}/?-*/{wildcards.ppr}/*{wildcards.ppr}.R {input.indir}/??-*/{wildcards.ppr}/*{wildcards.ppr}.R"
    
    echo "# Merged ceas r output between:" > {output.rcode}
    
    # Loop with ? and ?? regex to allow samples ordered correctly e.g. 1,2,...,10,11 and not 1,10,11,12,2...
    
    #
    # Header
    #
    for i in $SAMPLES
    do 
        head -17 $i >> {output.rcode}
    done
    
    echo "
    #
    # Average profile near TSS:
    #
    pdf('$OUTDIR/average_profile_near_tss_{wildcards.ppr}.pdf')
    
    x<-c(-3000.000000,-2950.000000,-2900.000000,-2850.000000,-2800.000000,-2750.000000,-2700.000000,-2650.000000,-2600.000000,-2550.000000,-2500.000000,-2450.000000,-2400.000000,-2350.000000,-2300.000000,-2250.000000,-2200.000000,-2150.000000,-2100.000000,-2050.000000,-2000.000000,-1950.000000,-1900.000000,-1850.000000,-1800.000000,-1750.000000,-1700.000000,-1650.000000,-1600.000000,-1550.000000,-1500.000000,-1450.000000,-1400.000000,-1350.000000,-1300.000000,-1250.000000,-1200.000000,-1150.000000,-1100.000000,-1050.000000,-1000.000000,-950.000000,-900.000000,-850.000000,-800.000000,-750.000000,-700.000000,-650.000000,-600.000000,-550.000000,-500.000000,-450.000000,-400.000000,-350.000000,-300.000000,-250.000000,-200.000000,-150.000000,-100.000000,-50.000000,0.000000,50.000000,100.000000,150.000000,200.000000,250.000000,300.000000,350.000000,400.000000,450.000000,500.000000,550.000000,600.000000,650.000000,700.000000,750.000000,800.000000,850.000000,900.000000,950.000000,1000.000000,1050.000000,1100.000000,1150.000000,1200.000000,1250.000000,1300.000000,1350.000000,1400.000000,1450.000000,1500.000000,1550.000000,1600.000000,1650.000000,1700.000000,1750.000000,1800.000000,1850.000000,1900.000000,1950.000000,2000.000000,2050.000000,2100.000000,2150.000000,2200.000000,2250.000000,2300.000000,2350.000000,2400.000000,2450.000000,2500.000000,2550.000000,2600.000000,2650.000000,2700.000000,2750.000000,2800.000000,2850.000000,2900.000000,2950.000000,3000.000000)
    
    y_merged <- c(
    " >> {output.rcode}
    
    # Get y values for each sample
          # Selecting line 342 and removing the "y<-c(" in the begining and susbstituting the ")" at the end of line with "," to create an R vector with y values for every samples.
    for i in $SAMPLES
    do
        head -342 $i | tail -1 | sed "s/y<-c(//" | sed "s/)$/,/">> {output.rcode}
    done
    
    # Substituting the last "," with ")" to close the R vector
    tail -1 {output.rcode} | sed "s/,$/)/" > {output.rcode}.tmp
    head -n-1 {output.rcode} > {output.rcode}.tmp2   # Remove the last line
    cat {output.rcode}.tmp2 {output.rcode}.tmp > {output.rcode}
    
    echo '
    library(RColorBrewer)
    colorPalette <- brewer.pal(12,"Paired")
    
    sample_names <- c("1-MNS-TGC-WT","2-MNS-TGC-KO","3-MNS-P-WT","4-MNS-P-KO","5-MNS-R-WT","6-MNS-R-KO","7-MNS-SC-WT","8-MNS-SC-KO","9-PSK-SC-WT","10-PSK-SC-KO","11-SNS-SC-WT","12-SNS-SC-KO")
    
    plot(x=c(-3000,3000), y=c(min(y_merged),max(y_merged)), type="n", main="Average Profile near TSS",xlab="Relative Distance to TSS (bp)",ylab="Average Profile")
    ' >> {output.rcode}
    
    COLOR=1
    for i in $SAMPLES; do
    
    echo '
    lines(x=x,y=' >> {output.rcode}
    
    head -342 $i | tail -1 | sed "s/y<-c(/c(/" >> {output.rcode}
    
    echo "
    , col=colorPalette[$COLOR])
    " >> {output.rcode}
    COLOR=$((COLOR+1))
    done
    
    echo '
    legend("topright", legend=sample_names, col=colorPalette, lty=1, lwd=2)
    abline(v=0.000000,lty=2,col=c("black"))
    dev.off()
    ' >> {output.rcode}
    
    
    echo "
    #
    # Average profile near TTS:
    #
    pdf('$OUTDIR/average_profile_near_tts_{wildcards.ppr}.pdf')
    
    y_merged <- c(
    " >> {output.rcode}
    
    # Get y values for each sample
          # Selecting line 346 and removing the "y<-c(" in the begining and susbstituting the ")" at the end of line with "," to create an R vector with y values for every samples.
    for i in $SAMPLES
    do
        head -346 $i | tail -1 | sed "s/y<-c(//" | sed "s/)$/,/">> {output.rcode}
    done
    
    # Substituting the last "," with ")" to close the R vector
    tail -1 {output.rcode} | sed "s/,$/)/" > {output.rcode}.tmp
    head -n-1 {output.rcode} > {output.rcode}.tmp2   # Remove the last line
    cat {output.rcode}.tmp2 {output.rcode}.tmp > {output.rcode}
    
    echo '
    plot(x=c(-3000,3000), y=c(min(y_merged),max(y_merged)), type="n", main="Average Profile near TTS",xlab="Relative Distance to TTS (bp)",ylab="Average Profile")
    ' >> {output.rcode}
    
    COLOR=1
    for i in $SAMPLES; do
    
    echo '
    lines(x=x,y=' >> {output.rcode}
    
    head -346 $i | tail -1 | sed "s/y<-c(/c(/" >> {output.rcode}
    
    echo "
    , col=colorPalette[$COLOR])
    " >> {output.rcode}
    COLOR=$((COLOR+1))
    done
    
    echo '
    legend("topright", legend=sample_names, col=colorPalette, lty=1, lwd=2)
    abline(v=0.000000,lty=2,col=c("black"))
    dev.off()
    ' >> {output.rcode}
    
    
    echo "
    #
    # Average Gene profile:
    #
    pdf('$OUTDIR/average_gene_profile_{wildcards.ppr}.pdf')
    
    x<-c(-1000.000000,-950.000000,-900.000000,-850.000000,-800.000000,-750.000000,-700.000000,-650.000000,-600.000000,-550.000000,-500.000000,-450.000000,-400.000000,-350.000000,-300.000000,-250.000000,-200.000000,-150.000000,-100.000000,-50.000000,0.000000,50.000000,100.000000,150.000000,200.000000,250.000000,300.000000,350.000000,400.000000,450.000000,500.000000,550.000000,600.000000,650.000000,700.000000,750.000000,800.000000,850.000000,900.000000,950.000000,1000.000000,1050.000000,1100.000000,1150.000000,1200.000000,1250.000000,1300.000000,1350.000000,1400.000000,1450.000000,1500.000000,1550.000000,1600.000000,1650.000000,1700.000000,1750.000000,1800.000000,1850.000000,1900.000000,1950.000000,2000.000000,2050.000000,2100.000000,2150.000000,2200.000000,2250.000000,2300.000000,2350.000000,2400.000000,2450.000000,2500.000000,2550.000000,2600.000000,2650.000000,2700.000000,2750.000000,2800.000000,2850.000000,2900.000000,2950.000000,3000.000000,3050.000000,3100.000000,3150.000000,3200.000000,3250.000000,3300.000000,3350.000000,3400.000000,3450.000000,3500.000000,3550.000000,3600.000000,3650.000000,3700.000000,3750.000000,3800.000000,3850.000000,3900.000000,3950.000000,4000.000000)
    
    y_merged <- c(
    " >> {output.rcode}
    
    # Get y values for each sample
          # Selecting line 350 and removing the "y<-c(" in the begining and susbstituting the ")" at the end of line with "," to create an R vector with y values for every samples.
    for i in $SAMPLES
    do
        head -350 $i | tail -1 | sed "s/y<-c(//" | sed "s/)$/,/">> {output.rcode}
    done
    
    # Substituting the last "," with ")" to close the R vector
    tail -1 {output.rcode} | sed "s/,$/)/" > {output.rcode}.tmp
    head -n-1 {output.rcode} > {output.rcode}.tmp2   # Remove the last line
    cat {output.rcode}.tmp2 {output.rcode}.tmp > {output.rcode}
    
    echo '
    
    plot(x=c(-1000,4000), y=c(min(y_merged),max(y_merged)), type="n", main="Average Gene profile", xlab="Upstream (bp), 3000 bp of Meta-gene, Downstream (bp)",ylab="Average Profile")
    ' >> {output.rcode}
    
    COLOR=1
    for i in $SAMPLES; do
    
    echo '
    lines(x=x,y=' >> {output.rcode}
    
    head -350 $i | tail -1 | sed "s/y<-c(/c(/" >> {output.rcode}
    
    echo "
    , col=colorPalette[$COLOR])
    " >> {output.rcode}
    COLOR=$((COLOR+1))
    done
    
    echo '
    legend("topright", legend=sample_names, col=colorPalette, lty=1, lwd=2)
    abline(v=0.000000,lty=2,col=c("black"))
    abline(v=3000.000000,lty=2,col=c("black"))
    dev.off()
    ' >> {output.rcode}
    
    
    echo "
    #
    # Average concatenated exon profile:
    #
    pdf('$OUTDIR/average_concatenated_exon_profile_{wildcards.ppr}.pdf')
    x<-c(0.000000,3.333333,6.666667,10.000000,13.333333,16.666667,20.000000,23.333333,26.666667,30.000000,33.333333,36.666667,40.000000,43.333333,46.666667,50.000000,53.333333,56.666667,60.000000,63.333333,66.666667,70.000000,73.333333,76.666667,80.000000,83.333333,86.666667,90.000000,93.333333,96.666667,100.000000)
    
    y_merged <- c(
    " >> {output.rcode}
    
    # Get y values for each sample
          # Selecting line 355 removing the "y<-c(" in the begining and susbstituting the ")" at the end of line with "," to create an R vector with y values for every samples.
    for i in $SAMPLES
    do
        head -355 $i | tail -1 | sed "s/y<-c(//" | sed "s/)$/,/">> {output.rcode}
    done
    
    # Substituting the last "," with ")" to close the R vector
    tail -1 {output.rcode} | sed "s/,$/)/" > {output.rcode}.tmp
    head -n-1 {output.rcode} > {output.rcode}.tmp2   # Remove the last line
    cat {output.rcode}.tmp2 {output.rcode}.tmp > {output.rcode}
    
    echo '
    plot(x=c(0,100), y=c(min(y_merged),max(y_merged)), type="n", main="Average Concatenated Exon Profile",xlab="Relative Location",ylab="Average Profile")
    ' >> {output.rcode}
    
    COLOR=1
    for i in $SAMPLES; do
    
    echo '
    lines(x=x,y=' >> {output.rcode}
    
    head -355 $i | tail -1 | sed "s/y<-c(/c(/" >> {output.rcode}
    
    echo "
    , col=colorPalette[$COLOR])
    " >> {output.rcode}
    COLOR=$((COLOR+1))
    done
    
    echo '
    legend("topright", legend=sample_names, col=colorPalette, lty=1, lwd=2)
    dev.off()
    ' >> {output.rcode}
    
    
    echo "
    #
    # Average concatenated intron profile:
    #
    pdf('$OUTDIR/average_concatenated_intron_profile_{wildcards.ppr}.pdf')
    
    y_merged <- c(
    " >> {output.rcode}
    
    # Get y values for each sample
          # Selecting line 358 removing the "y<-c(" in the begining and susbstituting the ")" at the end of line with "," to create an R vector with y values for every samples.
    for i in $SAMPLES
    do
        head -358 $i | tail -1 | sed "s/y<-c(//" | sed "s/)$/,/">> {output.rcode}
    done
    
    # Substituting the last "," with ")" to close the R vector
    tail -1 {output.rcode} | sed "s/,$/)/" > {output.rcode}.tmp
    head -n-1 {output.rcode} > {output.rcode}.tmp2   # Remove the last line
    cat {output.rcode}.tmp2 {output.rcode}.tmp > {output.rcode}
    
    echo '
    plot(x=c(0,100), y=c(min(y_merged),max(y_merged)), type="n", main="Average Concatenated Exon Profile",xlab="Relative Location",ylab="Average Profile")
    ' >> {output.rcode}
    
    COLOR=1
    for i in $SAMPLES; do
    
    echo '
    lines(x=x,y=' >> {output.rcode}
    
    head -358 $i | tail -1 | sed "s/y<-c(/c(/" >> {output.rcode}
    
    echo "
    , col=colorPalette[$COLOR])
    " >> {output.rcode}
    COLOR=$((COLOR+1))
    done
    
    echo '
    legend("topright", legend=sample_names, col=colorPalette, lty=1, lwd=2)
    dev.off()
    ' >> {output.rcode}
    
    rm -f $OUTDIR/*.tmp $OUTDIR/*.tmp2
    
    Rscript {output.rcode}
    """

rule zoom_on_merge_ceas_r_output:
    """
    The profiles are the same for positions, peaks and regions as the same wig is used for each.
    
    input commented to work on this rule on pedagogix. Remove comment to restore workflow functionality.
    """
    input:
        rcode="result/merge_ceas_r_output/{ppr}_code.r"
        # rcode=expand("result/ceas/danpos/{sample}/{ppr_in}/{sample}_{ppr_in}.R", ppr_in=PPR, sample=SAMPLES_MNASE)
    output:
        rcode="result/zoom_on_merge_ceas_r_output/{ppr}_code.r"
    shell:"""
    cat {input.rcode} | sed -e "s/y=c(min(y_merged),max(y_merged))/y=c(0,30)/g" \
        -e "s/merge_ceas_r_output/zoom_on_merge_ceas_r_output/g" > {output.rcode}
    
    Rscript {output.rcode}
    """

