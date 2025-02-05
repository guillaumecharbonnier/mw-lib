"""
Created:
    2017-03-09 10:35:41
Aim:
    This file will contain code using R to produce home-made analysis.
    out/ln/alias/experiments/MNase-Seq_mm10_peak_anno_stats/C_small_str_H2AL2_KO.txt out/ln/alias/experiments/MNase-Seq_mm10_peak_anno_stats/C_small_str.txt
"""

rule r_rscript:
    """
    Created:
        2018-12-03 15:28:55
        http://revigo.irb.hr/toR.jsp?table=1
    Test:
        out/ln/srf_from_src/r/script/revigo.R
        out/r/rscript/ln/updir/mw/src/r/script/test/done
    """
    input:
        script = 'out/{filler}.R',
        deps = r_smi_dep
    output:
        touch('out/{tool}{extra}/{filler}/done')
    wildcard_constraints:
        tool="r/rscript"
    params:
        extra = params_extra
    conda:
        "../envs/r_dev.yaml"
    shell:
        "Rscript {input.script} {params.extra}"

rule r_chen_sup_to_broad_or_sharp_tables:
    input:
        csv='out/bioconvert/xls2csv/wget/https/media.nature.com/original/nature-assets/ng/journal/v47/n10/extref/ng.3385-S4.csv',
        script='src/r/script/chen_sup_to_broad_or_sharp_tables.R'
    output:
        'out/r/chen_sup_to_broad_or_sharp_tables/broad.yaml'
    conda:
        "../envs/r_dev.yaml"
    shell:
        """
        Rscript {input.script}
        """

rule r_chen_test_clusterprofiler:
    input:
        yaml='out/r/chen_sup_to_broad_or_sharp_tables/broad.yaml',
        script='src/r/script/chen_test_clusterprofiler.R'
    output:
#        'out/r/chen_sup_to_broad_or_sharp_tables/broad.yaml'
    conda:
        "../envs/r_dev.yaml"
    threads:
        3 # One for each GO ontology
    shell:
        """
        Rscript {input.script}
        """

rule r_revigo_test:
    input:
        script="src/r/script/revigo.R"
    output:
        "out/r/revigo_test/scatterplot.pdf"
    conda:
        "../envs/r_dev.yaml"
    shell:
        """
        cp src/r/script/revigo.R out/r/revigo_test/revigo.R
        cd out/r/revigo_test
        Rscript {WDIR}/{input.script}
        mv Rplots.pdf scatterplot.pdf
        """

rule r_mnase_peak_anno_stats_to_mutiple_samples_barplots:
    """
    """
    input:
        peak_anno_stats = expand("out/ln/alias/experiments/MNase-Seq_mm10_peak_anno_stats/{sample}.txt", sample=[
            "naked_DNA_5min_r1",
            "naked_DNA_5min_r2",
            "naked_DNA_5min_r3",
            "naked_DNA_10min_r1",
            "naked_DNA_10min_r2",
            "naked_DNA_10min_r3",
            "S_b1_nuc",
            "S_b3_nuc",
            "S_b4_int_str",
            "S_b2_small_str",
            "S_b5_small_str",
            "C_small_str",
            "C_nuc",
            "C_small_str_H2AL2_KO",
            "C_nuc_H2AL2_KO",
            "P_nuc",
            "P_nuc_H2AL2_KO",
            "R_nuc",
            "R_nuc_H2AL2_KO"]),
        script="src/r/script/mnase_peak_anno_stats_to_mutiple_samples_barplots.R"
    output:
        "out/r/mnase_peak_anno_stats_to_mutiple_samples_barplots/log2_ratio_enrich_vs_uniform_distrib.pdf"
    conda:
        "../envs/r_chromstar.yaml"
    shell:
        """
        Rscript {input.script}
        """

rule r_multi_peak_anno_stats_to_heatmap:
    """
    Created:
        2018-11-14 22:08:22
    Aim:
    Test:
        out/r/multi_peak_anno_stats_to_heatmap/test19/log2_ratio_enrich_vs_uniform_distrib.pdf
    """
    input:
        peak_anno_stats_list = lambda wildcards: eval(mwconf['ids'][wildcards.peak_anno_stats_list_id]),
        script="src/r/script/multi_peak_anno_stats_to_heatmap.R"
    output:
        "out/r/multi_peak_anno_stats_to_heatmap/{peak_anno_stats_list_id}/log2_ratio_enrich_vs_uniform_distrib.pdf"
    conda:
        "../envs/r_greatr.yaml"
    shell:
        """
        FILES=`echo {input.peak_anno_stats_list} | sed 's/ /,/g'`
        OUTDIR=`dirname {output}`
        Rscript {input.script} -f $FILES -o $OUTDIR
        """



rule r_rmarkdown:
    """
    Created:
        2018-07-24 16:55:58
    Aim:
        Compile Rmd files with rmarkdown library

    out/r/rmarkdown/consensus_mhsc_t-cells_rgreat_revigo.html
    out/r/rmarkdown/cdr_hematopoiesis_rgreat_revigo.html
    """
    input:
        Rmd="src/rmd/{filler}.Rmd",
        script="src/r/script/rmarkdown.R"
    output:
        Rmd  = "out/r/rmarkdown/{filler}.Rmd",
        html = "out/r/rmarkdown/{filler}.html"#,
        #files= dynamic("out/r/rmarkdown/{{filler}}_files/figure-html/{doc}")
    conda:
        "../envs/r_rmd.yaml"
    shell:
        """
        # cp because I do not want to flood src directory with cache files.
        cp {input.Rmd} {output.Rmd}
        Rscript {input.script} -r {output.Rmd}
        """

rule r_rmarkdown_dynamic_enhancers_in_thymopoiesis:
    """
    Aim:
        Compile Rmd files with rmarkdown library

    out/r/rmarkdown/consensus_mhsc_t-cells_rgreat_revigo.html
    """
    input:
        Rmd="src/rmd/dynamic_enhancers_in_thymopoiesis.Rmd",
        script="src/r/script/rmarkdown.R",
        h3k27ac = "out/gtftk/coverage_--pseudo-count_0_--matrix-out_chrominfo-hg38-main-chr_bed-hg38-bs-hypometh-thymus-union-all_bw-hg38-H3K27ac-thymus-merged-wiq.txt",
        h3k27ac_ssw = "out/gtftk/coverage_--pseudo-count_0_--matrix-out_--upstream_2000_--downstream_2000_--nb-window_20_--n-highest_5_chrominfo-hg38-main-chr_bed-hg38-bs-hypometh-thymus-union-all_bw-hg38-H3K27ac-thymus-merged-wiq.txt",
        rna = 'out/gtftk/coverage_--upstream_10000_--downstream_10000_--pseudo-count_0_--matrix-out_chrominfo-hg38-main-chr_bed-hg38-bs-hypometh-thymus-union-all_bw-hg38-RNA-thymus.txt',
        rna_counts = "out/subread/featureCounts_-O_-t_exon_-g_gene_id_gtf-hg38-ensembl_bam-hg38-RNA-thymus.tsv",
        rna_counts_nodup = 'out/subread/featureCounts_-O_-t_exon_-g_gene_id_gtf-hg38-ensembl_bam-hg38-RNA-thymus-nodup.tsv',
        rmsk = "out/awk/extract_bed6_from_rmsk/gunzip/to-stdout/wget/http/hgdownload.cse.ucsc.edu/goldenpath/hg38/database/rmsk.bed",
        dist_tss = "out/bedtools/closest_-d_-t_first_bed-hg38-ensembl-r93-tss-protein-coding-TR-IG/sort/coordinates-bed/bedtools/merge/bedtools/multiinter_bed-hg38-bs-hypometh-thymus/multiinter.bed",
        atac = 'out/gtftk/coverage_--pseudo-count_0_--matrix-out_chrominfo-hg38-main-chr_bed-hg38-bs-hypometh-thymus-union-all_bw-hg38-ATAC-thymus-wiq.txt',
        closest_upstream = 'out/bedtools/closest_-D_ref_-id_-k_1_-t_first_bed-hg38-ensembl-r93-tss-protein-coding-TR-IG/sort/coordinates-bed/bedtools/merge/bedtools/multiinter_bed-hg38-bs-hypometh-thymus/multiinter.bed',
        closest_downstream = 'out/bedtools/closest_-D_ref_-iu_-k_1_-t_first_bed-hg38-ensembl-r93-tss-protein-coding-TR-IG/sort/coordinates-bed/bedtools/merge/bedtools/multiinter_bed-hg38-bs-hypometh-thymus/multiinter.bed',
        features = "out/bedtools/merge/bedtools/multiinter_bed-hg38-bs-hypometh-thymus/multiinter.bed",
        bed_paths = ["out/awk/extract_xls_coordinates_to_bed6/macs2/callpeak_--gsize_hs_--nomodel/ln/alias/experiments/hg38_H3K27ac_thymus/CD34_over_input_peaks.bed",
        "out/awk/extract_xls_coordinates_to_bed6/macs2/callpeak_--gsize_hs_--nomodel/ln/alias/experiments/hg38_H3K27ac_thymus/EC_over_input_peaks.bed",
        "out/awk/extract_xls_coordinates_to_bed6/macs2/callpeak_--gsize_hs_--nomodel/ln/alias/experiments/hg38_H3K27ac_thymus/LC_over_input_peaks.bed",
        "out/awk/extract_xls_coordinates_to_bed6/macs2/callpeak_--gsize_hs_--nomodel/ln/alias/experiments/hg38_H3K27ac_thymus/SP4_over_input_peaks.bed",
        "out/awk/extract_xls_coordinates_to_bed6/macs2/callpeak_--gsize_hs_--nomodel/ln/alias/experiments/hg38_H3K27ac_thymus/SP8_over_input_peaks.bed"],
        atac_paths = ["out/awk/extract_xls_coordinates_to_bed6/macs2/noctrl_callpeak_--broad_--gsize_hs_--keep-dup_all/ln/alias/experiments/hg38_ATAC_thymus/CD34_peaks.bed",
        "out/awk/extract_xls_coordinates_to_bed6/macs2/noctrl_callpeak_--broad_--gsize_hs_--keep-dup_all/ln/alias/experiments/hg38_ATAC_thymus/EC_peaks.bed",
        "out/awk/extract_xls_coordinates_to_bed6/macs2/noctrl_callpeak_--broad_--gsize_hs_--keep-dup_all/ln/alias/experiments/hg38_ATAC_thymus/LC_peaks.bed"],
        cpg_paths = ['out/ln/alias/experiments/BS_Blueprint_healthy_methylation_calls/TH_CD34_110.bw',
        'out/ln/alias/experiments/BS_Blueprint_healthy_methylation_calls/TH_EC_91.bw',
        'out/ln/alias/experiments/BS_Blueprint_healthy_methylation_calls/TH_EC_101.bw',
        'out/ln/alias/experiments/BS_Blueprint_healthy_methylation_calls/TH_LC_91.bw',
        'out/ln/alias/experiments/BS_Blueprint_healthy_methylation_calls/TH_LC_101.bw',
        'out/ln/alias/experiments/BS_Blueprint_healthy_methylation_calls/TH_CD4_91.bw',
        'out/ln/alias/experiments/BS_Blueprint_healthy_methylation_calls/TH_CD4_101.bw',
        'out/ln/alias/experiments/BS_Blueprint_healthy_methylation_calls/TH_CD8_91.bw',
        'out/ln/alias/experiments/BS_Blueprint_healthy_methylation_calls/TH_CD8_101.bw'],
        intersect_merge_and_multiinter = 'out/bedtools/intersect_bs_hypometh_thymus_merge_and_multiinter/test.bed'
    output:
        Rmd  = "out/r/dynamic_enhancers_in_thymopoiesis.Rmd",
        html = "out/r/dynamic_enhancers_in_thymopoiesis.html"#,
        #files= dynamic("out/r/rmarkdown/{{filler}}_files/figure-html/{doc}")
    conda:
        "../envs/r_rmd.yaml"
    shell:
        """
        # cp because I do not want to flood src directory with cache files.
        cp {input.Rmd} {output.Rmd}
        Rscript {input.script} -r {output.Rmd}
        """

# Dynamic is removed from snakemake v8.
# rule r_rmarkdown_dynamic_consensus_mhsc_tcells_rgreat_revigo:
#     """
#     Created:
#         2018-07-24 16:55:58
#     Aim:
#         Compile Rmd files with rmarkdown library.
#         Dynamic allow to track all files produced by knitr for this Rmd.
#     Test:
#         out/r/rmarkdown_dynamic/consensus_mhsc_t-cells_rgreat_revigo.html
#     """
#     input:
#         Rmd="src/rmd/consensus_mhsc_t-cells_rgreat_revigo.Rmd",
#         #Rmd="src/rmd/{filler}.Rmd",
#         Rscript="opt/miniconda/envs/r/bin/Rscript",
#         #Rscript="opt/miniconda/envs/rgreat_revigo/bin/Rscript",
#         script="src/r/script/rmarkdown.R"
#     output:
#         files= dynamic("out/r/rmarkdown_dynamic/consensus_mhsc_t-cells_rgreat_revigo{files}")
#         #files= dynamic("out/r/rmarkdown_dynamic/consensus_mhsc_tcells_rgreat_revigo_files/figure-html/{doc}")
#     params:
#         Rmd  = "out/r/rmarkdown_dynamic/consensus_mhsc_t-cells_rgreat_revigo.Rmd",
#         html = "out/r/rmarkdown_dynamic/consensus_mhsc_t-cells_rgreat_revigo.html",
#     shell:
#         """
#         set +u; source opt/miniconda/bin/activate r; set -u
#         # cp because I do not want to flood src directory with cache files.
#         cp {input.Rmd} {params.Rmd}
#         ln -sf ../../../out out/r/rmarkdown_dynamic/out
#         Rscript {input.script} -r {params.Rmd}
#         #{input.Rscript} {input.script} -r {params.Rmd}
#         """

rule r_bookdown:
    """
    Created:
        2018-07-24 18:13:50
    Aim:
        Compile Rmd files with bookdown library
    Test:
        out/r/bookdown/bookdown-demo-master/index.html
    """
    input:
        Rmd="src/rmd/{filler}/index.Rmd",
        Rscript="opt/miniconda/envs/r/bin/Rscript",
        script="src/r/script/bookdown.R"
    output:
        Rmd  = "out/r/bookdown/{filler}/index.Rmd",
        html = "out/r/bookdown/{filler}/index.html"
    shell:
        """
        set +u; source opt/miniconda/bin/activate r; set -u
        # cp because I do not want to flood src directory with cache files.
        cp -r $(dirname {input.Rmd}) $(dirname $(dirname {output.Rmd}))
        cd `dirname {output.Rmd}`
        {WDIR}/{input.Rscript} {WDIR}/{input.script} -r index.Rmd
        """

rule r_chromstar:
    """
    Created:
        2018-10-21 16:20:20
    Aim:
    """
    input:
        "out/wget/https/bioconductor.org/packages/release/bioc/vignettes/chromstaR/inst/doc/chromstaR.R"

rule r_fragment_distrib_from_samtools_view_get_tlen:
    """
    Created:
        2018-07-02 11:49:59
    """
    input:
        txt=expand("out/samtools/view_get_tlen/ln/alias/experiments/mm10_H4K5ac_comparison/{treatment}-{condition}_{run}.txt", treatment=["H4K5ac","Input"], condition=["Nut-WT","Nut-KO"], run=["run207","run215", "run140-141"]),
        code="src/r/script/fragment_distrib_from_samtools_view_get_tlen.R"
    shell:
        """

        """

rule r_fragment_distrib_from_samtools_view_get_tlen_fastkd1:
    """
    Created:
        2018-08-23 11:50:40
    Aim:

    """
    input:
        txt=expand("out/samtools/view_get_tlen/ln/alias/experiments/ChIP-Seq_Spike-in_REH/GRCh38/run{run}-{chip}-{cellLine}.txt",run=["219","223"], chip=["H4K5ac","H4K5cr","Input"], cellLine=["gev","g3","g5"]),
        txt_bu_219=expand("out/samtools/view_get_tlen/ln/alias/experiments/ChIP-Seq_Spike-in_REH/GRCh38/run219-H4K5bu-{cellLine}.txt", cellLine=["gev","g3","g5f"]),
        txt_bu_223=expand("out/samtools/view_get_tlen/ln/alias/experiments/ChIP-Seq_Spike-in_REH/GRCh38/run223-H4K5bu-{cellLine}.txt", cellLine=["gev","g3","g5f","g5"]),
        code="src/r/script/fragment_distrib_from_samtools_view_get_tlen_fastkd1.R",
        Rscript="opt/miniconda/envs/r/bin/Rscript"
    shell:
        """
        {input.Rscript} {input.code}
        """

rule r_great_heatmap_from_first_clustering:
    """
    Created:
        2018-10-01 13:40:19
    Test:
        out/r/great_heatmap_bed-hg19-cpg-meth-call-clusters/enrichment_tables.Rdata
        out/r/great_heatmap_g-GOBP_r-NULL_bed-hg19-h3k27ac-on-low-cpg-meth-call-clusters/enrichment_tables.Rdata
        out/r/great_heatmap__g-GOBP_m-bed-hg19-alex-test/enrichment_tables.done
        out/r/great_heatmap__g-GOBP__m-Hyper_Fold_Enrichment,Binom_Bonf_PValue,Hyper_Bonf_PValue,Post_Filter_Binom_Rank__s-greater,lower,lower,lower__t-2,0.05,0.05,5_bed-hg19-alex-test/enrichment_tables.done
        out/r/great_heatmap__g-MSDBP__m-Binom_Fold_Enrichment,Binom_Bonf_PValue,Hyper_Bonf_PValue__s-greater,lower,lower__t-1.5,0.001,0.001_bed-hg19-cpg-meth-call-clusters/enrichment_tables.done
    """
    input:
        bed_list = lambda wildcards: eval(mwconf['ids'][wildcards.bed_list_id]),
        script="src/r/pkg/gcfun/R/GREATR.R",
        #Rscript="opt/miniconda/envs/r_great_heatmap/bin/Rscript"
    output:
        #tables="out/r/great_heatmap{extra}_bed-{bed_list_id}/enrichment_tables.Rdata"
        done="out/r/great_heatmap{extra}_bed-{bed_list_id}/enrichment_tables.done",
        pdf=expand("out/r/great_heatmap{{extra}}_bed-{{bed_list_id}}/{plot}.pdf", plot=["facet_fold_enrich"])
    params:
        #indir="out/crossmap/hg38_to_hg19/r/dynamic_enchancers_in_thymopoiesis/dClust/rowFeature-no_rmsk_mxy__no_donor_effect__distal/sortingMethod-kmeans_centers-8_nstart-100_itermax-200000_algorithm-Lloyd/sortingFeature-cpg_meth_call/subSortingMethod-kmeans_centers-8_nstart-100_itermax-200000_algorithm-Lloyd/subSortingFeature-H3K27ac_peaks.ATAC_peaks",
        #outdir="out/r/great_heatmap/crossmap/hg38_to_hg19/r/dynamic_enchancers_in_thymopoiesis/dClust/rowFeature-no_rmsk_mxy__no_donor_effect__distal/sortingMethod-kmeans_centers-8_nstart-100_itermax-200000_algorithm-Lloyd/sortingFeature-cpg_meth_call/subSortingMethod-kmeans_centers-8_nstart-100_itermax-200000_algorithm-Lloyd/subSortingFeature-H3K27ac_peaks.ATAC_peaks"
    conda:
        "../envs/r_great_heatmap.yaml"
    shell:
        """
        EXTRA=`echo {wildcards.extra} | sed -e 's/-/ /g' -e 's/__/ -/g'`
        FILES=`echo {input.beds} | sed 's/ /,/g'`
        OUTDIR=`dirname {output.done}`
        echo "OUTDIR: $OUTDIR"
        Rscript {input.script} -f $FILES -o $OUTDIR $EXTRA

        touch {output.done}
        """

rule r_fraction_of_genome_above_n:
    input:
        bw=expand("out/deepTools/bamCoverage_binSize-{binSize}_normalizeUsing-RPKM_extendReads-/ln/alias/experiments/mm10_H4K5ac_comparison/{treatment}-{condition}_{run}.bw", treatment=["H4K5ac","Input"], condition=["Nut-WT","Nut-KO"], run=["run207","run215", "run140-141"], binSize=["200","1000","10000"]),



#def r_data_dep(wildcards):
#    """
#    Created:
#        2017-03-24 14:35:49
#    Aim:
#        It could be nice to have as input of the latex rule the list of files needed in the tex file to compile.
#    """
#    id=wildcards['id']
#    rcode="src/r/script/"+id+".R"
#
#    # Meaning of the regex:
#    # \s*[^%] means line that is not commented in latex.
#    # \\includegraphics.+{ means includegraphic call with options until start of file path.
#    # (.*) is the file path extracted.
#    # } is the end of the includegraphic call.
#    # This regex implies avoiding additionnal {} in the file path.
#    pattern = re.compile(r'read\.table\(file="(.*)"')
#    paths = []
#
#    with open (rcode, "rt") as infile:
#        for line in infile:
#            print(line)
#            m = pattern.search(line)
#            if m:
#                path = m.group(1)
#                paths.append(path)
#
#    return(paths)
#
#def r_data_out(wildcards):
#    """
#    Created:
#        2017-03-24 15:03:18
#    Note: Can not work because snakemake accepts only functions as inputs and not outputs.
#    Aim:
#        It could be nice to have as input of the latex rule the list of files needed in the tex file to compile.
#    """
#    id=wildcards['id']
#    rcode="src/r/script/"+id+".R"
#
#    # Meaning of the regex:
#    # \s*[^%] means line that is not commented in latex.
#    # \\includegraphics.+{ means includegraphic call with options until start of file path.
#    # (.*) is the file path extracted.
#    # } is the end of the includegraphic call.
#    # This regex implies avoiding additionnal {} in the file path.
#    pattern = re.compile(r'write\.table\(file="(.*)"')
#    paths = []
#
#    with open (rcode, "rt") as infile:
#        for line in infile:
#            print(line)
#            m = pattern.search(line)
#            if m:
#                path = m.group(1)
#                paths.append(path)
#
#    return(paths)
#
#rule r:
#    """
#    Created:
#        2017-03-24 12:09:24
#    Aim:
#        General rule that execute a Rscript and put the results in a specific folder.
#    Note:
#        This idea can not work because of this:
#        SyntaxError:
#        Only input files can be specified as functions
#    Test:
#        out/r/script/test.done
#    """
#    input:
#        rscript="opt/miniconda/envs/r/bin/Rscript",
#        code="src/r/script/{id}.R",
#        dep=r_data_dep
#    output:
#        done="out/r/script/{id}.done",
#        dep=r_data_out
#    shell:
#        """
#        {input.rscript} {input.code}
#        """

rule r_beanplot_on_gtftk_coverage_tss_brdt_brd4_h4k5ac_classes:
    """
    Created:
        2017-03-09 10:37:36
    Aim:
        Take gtftk coverage matrices as input to produces violin plots showing the distribution of acetylation on each class of TSS.
    Note:
        For the moment it is a draft for the Nut analysis and produces other distribution plots too.
    """
    input:
        rscript="opt/miniconda/envs/r/bin/Rscript",
        code="src/r/script/beanplot.R",
        ##First test:
        #txt_tss_with_brdt_and_h4k5ac="out/gtftk/coverage_bed_test_tss_with_brdt_and_h4k5ac.txt",
        #txt_tss_with_brd4_and_h4k5ac="out/gtftk/coverage_bed_test_tss_with_brd4_and_h4k5ac.txt",
        #txt_tss_without_brdt_brd4_and_with_h4k5ac="out/gtftk/coverage_bed_test_tss_without_brdt_brd4_and_with_h4k5ac.txt",
        #txt_tss_with_h4k5ac="out/gtftk/coverage_bed_test_tss_with_h4k5ac.txt",
        ##Second test with RPKM normalization:
        #txt_tss_with_brdt_and_h4k5ac="out/gtftk/coverage_bed_test_tss_with_brdt_and_h4k5ac_minMappingQuality0.txt",
        #txt_tss_with_brd4_and_h4k5ac="out/gtftk/coverage_bed_test_tss_with_brd4_and_h4k5ac_minMappingQuality0.txt",
        #txt_tss_without_brdt_brd4_and_with_h4k5ac="out/gtftk/coverage_bed_test_tss_without_brdt_brd4_and_with_h4k5ac_minMappingQuality0.txt",
        #txt_tss_with_h4k5ac="out/gtftk/coverage_bed_test_tss_with_h4k5ac_minMappingQuality0.txt",
        #txt_brdt_peaks="out/gtftk/coverage_bed_test_brdt_peaks_minMappingQuality0.txt",
        #txt_brd4_peaks="out/gtftk/coverage_bed_test_brd4_peaks_minMappingQuality0.txt",
        #txt_h4k5ac_peaks="out/gtftk/coverage_bed_test_h4k5ac_peaks_minMappingQuality0.txt"
        #Third test with SES:
        txt_tss_with_brdt_and_h4k5ac              = "out/gtftk/coverage/bedtss_with_brdt_and_h4k5ac_bwNut_H4K5ac_bamCompare_SES.txt",
        txt_tss_with_brd4_and_h4k5ac              = "out/gtftk/coverage/bedtss_with_brd4_and_h4k5ac_bwNut_H4K5ac_bamCompare_SES.txt",
        txt_tss_without_brdt_brd4_and_with_h4k5ac = "out/gtftk/coverage/bedtss_without_brdt_brd4_and_with_h4k5ac_bwNut_H4K5ac_bamCompare_SES.txt",
        txt_tss_with_h4k5ac                       = "out/gtftk/coverage/bedtss_with_h4k5ac_bwNut_H4K5ac_bamCompare_SES.txt",
        txt_brdt_peaks                            = "out/gtftk/coverage/bedbrdt_peaks_bwNut_H4K5ac_bamCompare_SES.txt",
        txt_brd4_peaks                            = "out/gtftk/coverage/bedbrd4_peaks_bwNut_H4K5ac_bamCompare_SES.txt",
        txt_h4k5ac_peaks                          = "out/gtftk/coverage/bedh4k5ac_peaks_bwNut_H4K5ac_bamCompare_SES.txt",
        txt_tss_nut_ko_downreg                    = "out/gtftk/coverage/bedtss_nut_ko_downreg_bwNut_H4K5ac_bamCompare_SES.txt",
        txt_tss_nut_ko_upreg                      = "out/gtftk/coverage/bedtss_nut_ko_upreg_bwNut_H4K5ac_bamCompare_SES.txt"
    output:
        pdf_raw="out/r/beanplot_on_gtftk_coverage_brdt_brd4_h4k5ac_classes_raw_ses.pdf",
        pdf_ratio="out/r/beanplot_on_gtftk_coverage_brdt_brd4_h4k5ac_classes_ratio_ses.pdf"
    shell:
        """
        {input.rscript} {input.code} \
            --input {input.txt_tss_with_brdt_and_h4k5ac} \
                    {input.txt_tss_with_brd4_and_h4k5ac} \
                    {input.txt_tss_without_brdt_brd4_and_with_h4k5ac} \
                    {input.txt_tss_with_h4k5ac} {input.txt_brdt_peaks} \
                    {input.txt_brd4_peaks} \
                    {input.txt_h4k5ac_peaks} \
                    {input.txt_tss_nut_ko_downreg} \
                    {input.txt_tss_nut_ko_upreg} \
            --output-raw {output.pdf_raw} \
            --output-ratio {output.pdf_ratio}
        """
    #run:
    #    R("""
    #    library(beanplot)

    #    txt_tss_with_brdt_and_h4k5ac <- read.table(
    #        file='{input.txt_tss_with_brdt_and_h4k5ac}',
    #        comment.char = "#",
    #        header = T)

    #    txt_tss_with_brd4_and_h4k5ac <- read.table(
    #        file='{input.txt_tss_with_brd4_and_h4k5ac}',
    #        comment.char = "#",
    #        header = T)

    #    txt_tss_without_brdt_brd4_and_with_h4k5ac <- read.table(
    #        file='{input.txt_tss_without_brdt_brd4_and_with_h4k5ac}',
    #        comment.char = "#",
    #        header = T)

    #    # Having trouble sometimes with default bw settings:
    #    # http://r.789695.n4.nabble.com/beanplot-Error-sample-is-too-sparse-to-find-TD-td4291739.html

    #    pdf('{output.pdf}')
    #    beanplot(
    #        txt_tss_with_brdt_and_h4k5ac$H4K5ac-Nut-WT,
    #        txt_tss_with_brdt_and_h4k5ac$H4K5ac-Nut-KO,
    #        txt_tss_with_brd4_and_h4k5ac$H4K5ac-Nut-WT,
    #        txt_tss_with_brd4_and_h4k5ac$H4K5ac-Nut-KO,
    #        txt_tss_without_brdt_brd4_and_with_h4k5ac$H4K5ac-Nut-WT,
    #        txt_tss_without_brdt_brd4_and_with_h4k5ac$H4K5ac-Nut-KO,
    #        bw="nrd0",
    #        side='both',
    #        names=c('+Brdt+H4K5ac','+Brd4+H4K5ac','-Brdt-Brd4+H4K5ac'),
    #        ylab='RPKM coverage in TSS',
    #        xlab='WT (half left) versus KO (half right)',
    #        #main=paste0('n peaks WT: ', dim(xls1)[1], '; KO: ', dim(xls2)[1], '.')
    #        )
    #    dev.off()
    #    """)

rule r_format_featureCounts_to_cluster_pseudoLog2Count:
    """
    Created:
        2017-03-28 16:48:00
    Aim:
    Note:
        out/r/format_featureCounts_to_cluster_pseudoLog2Count/featureCounts/gtfcuffmerge-GRCm38-H2AL2-outFilterMultimapNmax-1000_texon_ggene_id_h2al2_rnaseq.tsv
    """
    input:
        tsv="out/{filler}.tsv"
    output:
        tsv="out/r/format_featureCounts_to_cluster_pseudoLog2Count/{filler}.tsv"
    run:
        R("""
        data <- read.table(
            file='{input.tsv}',
            #file='out/featureCounts/gtfcuffmerge-GRCm38-H2AL2-outFilterMultimapNmax-1000_texon_ggene_id_h2al2_rnaseq.tsv',
            comment.char = "#",
            header = T)
        print(colnames(data))

        #data <- data[,-c('Chr','Start','End','Strand','Length')]
        data_counts <- data[,-c(1:6)]

        data_counts <- log2(data_counts + 1)

        data_out <- data.frame(data[,'Geneid'], data_counts)

        write.table(
            x = data_out,
            file = '{output.tsv}',
            #file = "test.tsv",
            append = FALSE,
            quote = FALSE,
            sep = '\t',
            row.names=FALSE)

        # Should be faster but require library(data.table)
        # http://stackoverflow.com/questions/7072159/how-do-you-remove-columns-from-a-data-frame
        #[,'Chr':=NULL]
        #data[,Start:=NULL]
        #data[,End:=NULL]
        #data[,Strand:=NULL]
        #data[,Length:=NULL]



        print(colnames(data))

        """)
