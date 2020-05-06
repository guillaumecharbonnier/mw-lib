# 2019-02-04 18:22:12
# The whole file is commented because a new peak_anno is likely to come bundled with pygtftk soon.
# and current version have multiple broken points (eg no more input_gtf)
#rule gtftk_peak_anno_chrominfo_gtf:
#    """
#    Created: 
#        2017-03-07 16:19:51
#    Modified:
#        2017-06-09 12:47:04 - changed 16 threads to 1 since memory usage is reasonable now.
#    Aim:
#        Test as an alternative to featureCounts. Work only for mm10 right now for tests.
#    Note:
#        Untested
#    Test:
#        out/gtftk/peak_anno_chrominfo-mm10-main-chr_gtf-GRCm38-main-chr/awk/extract_main_chr/sh/danpos_xls_to_bed3/danpos/dtriple_v4/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run113_run125_run126/PSK-SC-WT/danpos.smooth.positions/00_peak_anno_diagrams.pdf
#    """
#    input:
#        chrominfo=input_chrominfo,
#        gtf = lambda wildcards: eval(config['ids'][wildcards.gtf_id]),
#        bed="out/{filler}.bed"
#    output:
#        pdf       = "out/gtftk/peak_anno_chrominfo-{chrominfo_id}_{gtf_id}/{filler}/00_peak_anno_diagrams.pdf"
#    params:
#        outputdir = "out/gtftk/peak_anno_chrominfo-{chrominfo_id}_{gtf_id}/{filler}"
#    threads:
#        1
#    wildcard_constraints:
#        gtf_id="[a-zA-Z0-9-]+",
#        chrominfo_id="[a-zA-Z0-9-]+"
#    conda:
#        "../envs/pygtftk.yaml"
#    shell:
#        """
#        gtftk peak_anno \
#            --inputfile {input.gtf} \
#            --outputdir {params.outputdir} \
#            --chrom-info {input.chrominfo} \
#            --peak-file {input.bed} \
#            --no-date
#        """
#
#rule gtftk_peak_anno_chrominfo_gtf_bed_anno_list_more_bed_labels:
#    """
#    Created: 
#        2017-03-07 16:19:51
#    Modified:
#        2017-06-09 12:47:04 - changed 16 threads to 1 since memory usage is reasonable now.
#    Aim:
#        Test as an alternative to featureCounts. Work only for mm10 right now for tests.
#    Note:
#        It seems that R fails to produce correct plots if the file path is longer thant 260. A simple way to prevent that inside mw long file tree is to change directory to the output one before launching gtftk so hardcoded paths in R code will be shorts. Need to test this assumption.
#    Test:
#        out/gtftk/peak_anno_chrominfo-mm10-main-chr_gtf-mm10-main-chr_bed-anno-list-mm10-main-chr-gsat-mm-synrep-mm_more-bed-labels/inp/unit_tests/gtftk_peak_anno_chrominfo_gtf_bed_anno_list_more_bed_labels/peak-file1.bed
#        out/gtftk/peak_anno_chrominfo-mm10-main-chr_gtf-mm10-main-chr_bed-anno-list-mm10-main-chr-tss-classes-and-brd-peaks_more-bed-labels/awk/extract_main_chr/awk/convert_bedpe_to_bed6_insert_size/bedtools/bamtobed_bedpe/samtools/sort_n/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run140_run141/H4K5ac-Nut-KO/00_peak_anno_diagrams.pdf
#        out/gtftk/peak_anno_chrominfo-mm10-main-chr_gtf-GRCm38-main-chr_bed-anno-list-mm10-main-chr-gsat-mm-synrep-mm_more-bed-labels/awk/extract_main_chr/awk/convert_bedpe_to_bed6_insert_size/bedtools/bamtobed_bedpe/samtools/sort_n/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-30_lmax-100/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band2/00_peak_anno_diagrams.pdf
#        out/gtftk/peak_anno_chrominfo-mm10-main-chr_gtf-GRCm38-main-chr_bed-anno-list-mm10-main-chr-nuc-r-cpg-h2alap1-peaks_more-bed-labels/awk/extract_from_bed4_feature-E12/bedtools/makewindows_w-200_i-src/ChromHMM/LearnModel_test11_numstates-20_assembly-mm10/spermatozoa_20_segments/00_peak_anno_diagrams.pdf
#        2017-09-28 10:01:34:
#            out/gtftk/peak_anno_chrominfo-mm10-main-chr_gtf-mm10-main-chr_bed-anno-list-mm10-main-chr-all-2017-09-28_more-bed-labels/awk/extract_main_chr/awk/convert_bedpe_to_bed6_insert_size/bedtools/bamtobed_bedpe/samtools/sort_n/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run140_run141/H4K5ac-Nut-KO/00_peak_anno_diagrams.pdf
#
#    """
#    input:
#        gtftk="opt/miniconda/envs/gtftk_dev/bin/gtftk",
#        chrominfo=input_chrominfo,
#        gtf=input_gtf,
#        bed_anno_list=input_bed_anno_list_gtftk_peak_anno,
#        bed="out/{filler}.bed"
#    output:
#        pdf=      "out/gtftk/peak_anno_chrominfo-{chrominfo_id}_gtf-{gtf_id}_bed-anno-list-{bed_anno_list_id}_more-bed-labels/{filler}/00_peak_anno_diagrams.pdf",
#        stats=      "out/gtftk/peak_anno_chrominfo-{chrominfo_id}_gtf-{gtf_id}_bed-anno-list-{bed_anno_list_id}_more-bed-labels/{filler}/00_peak_anno_stats.txt"
#    params:
#        outputdir="out/gtftk/peak_anno_chrominfo-{chrominfo_id}_gtf-{gtf_id}_bed-anno-list-{bed_anno_list_id}_more-bed-labels/{filler}",
#        tmpdir=   "out/gtftk/peak_anno_chrominfo-{chrominfo_id}_gtf-{gtf_id}_bed-anno-list-{bed_anno_list_id}_more-bed-labels/{filler}/tmp",
#        more_bed_labels=params_more_bed_labels_gtftk_peak_anno
#    threads:
#        1
#    wildcard_constraints:
#        gtf_id="[a-zA-Z0-9-]+",
#        chrominfo_id="[a-zA-Z0-9-]+",
#        bed_anno_list_id="[a-zA-Z0-9-]+",
#        #more_bed_labels_id="[a-zA-Z0-9-]+",
#    shell:
#        """
#        # Needed to get bedtools in the PATH:
#        export PATH="{WDIR}/opt/miniconda/envs/gtftk_dev/bin/:$PATH"
#        
#        wc -l {input.bed_anno_list}
#
#        cd {params.outputdir}
#
#        #echo '{input.bed_anno_list}'
#
#        INPUT_BED_ANNO_LIST=`echo '{input.bed_anno_list}'| sed -e 's| | {WDIR}/|g' -e "s|^|{WDIR}/|g"`
#
#        #echo INPUT_BED_ANNO_LIST:
#        echo $INPUT_BED_ANNO_LIST
#
#        {WDIR}/{input.gtftk} peak_anno \
#            --inputfile {WDIR}/{input.gtf} \
#            --more-bed $INPUT_BED_ANNO_LIST \
#            --more-bed-labels {params.more_bed_labels} \
#            --outputdir . \
#            --chrom-info {WDIR}/{input.chrominfo} \
#            --peak-file {WDIR}/{input.bed} \
#            --no-date
#        """
#
#rule gtftk_peak_anno_chrominfo_gtf_bed_anno_list_more_bed_labels_pdf_width_pdf_height:
#    """
#    Created: 
#        2017-10-04 16:10:56
#    Aim:
#        Since a recent commit, gtftk peak_anno give free dimension for the plot by default which is bad if you want to produces article-ready plots for a long list of features.
#    Test:
#    """
#    input:
#        gtftk="opt/miniconda/envs/gtftk_dev/bin/gtftk",
#        chrominfo=input_chrominfo,
#        gtf=input_gtf,
#        bed_anno_list=input_bed_anno_list_gtftk_peak_anno,
#        bed="out/{filler}.bed"
#    output:
#        pdf=      "out/gtftk/peak_anno_chrominfo-{chrominfo_id}_gtf-{gtf_id}_bed-anno-list-{bed_anno_list_id}_more-bed-labels_pdf-width-{pdf_width}_pdf-height-{pdf_height}/{filler}/00_peak_anno_diagrams.pdf"
#    params:
#        outputdir="out/gtftk/peak_anno_chrominfo-{chrominfo_id}_gtf-{gtf_id}_bed-anno-list-{bed_anno_list_id}_more-bed-labels_pdf-width-{pdf_width}_pdf-height-{pdf_height}/{filler}",
#        tmpdir=   "out/gtftk/peak_anno_chrominfo-{chrominfo_id}_gtf-{gtf_id}_bed-anno-list-{bed_anno_list_id}_more-bed-labels_pdf-width-{pdf_width}_pdf-height-{pdf_height}/{filler}/tmp",
#        more_bed_labels=params_more_bed_labels_gtftk_peak_anno
#    threads:
#        1
#    wildcard_constraints:
#        gtf_id="[a-zA-Z0-9-]+",
#        chrominfo_id="[a-zA-Z0-9-]+",
#        bed_anno_list_id="[a-zA-Z0-9-]+",
#        pdf_width="[0-9]+",
#        pdf_height="[0-9]+"
#        #more_bed_labels_id="[a-zA-Z0-9-]+",
#    shell:
#        """
#        # Needed to get bedtools in the PATH:
#        export PATH="{WDIR}/opt/miniconda/envs/gtftk_dev/bin/:$PATH"
#
#        cd {params.outputdir}
#
#        INPUT_BED_ANNO_LIST=`echo '{input.bed_anno_list}'| sed -e 's| | {WDIR}/|g' -e "s|^|{WDIR}/|g"`
#
#        {WDIR}/{input.gtftk} peak_anno \
#            --inputfile {WDIR}/{input.gtf} \
#            --more-bed $INPUT_BED_ANNO_LIST \
#            --more-bed-labels '{params.more_bed_labels}' \
#            --outputdir . \
#            --chrom-info {WDIR}/{input.chrominfo} \
#            --peak-file {WDIR}/{input.bed} \
#            --pdf-width {wildcards.pdf_width} \
#            --pdf-height {wildcards.pdf_height} \
#            --no-date
#        """
#
#rule gtftk_peak_anno_chrominfo_gtf_more_keys:
#    """
#    Created: 
#        2017-06-26 15:47:51
#    Modified:
#    Aim:
#        Test --more-key feature.
#    Test:
#        out/gtftk/peak_anno_chrominfo-mm10-main-chr_gtf-GRCm38-main-chr_more-keys-gene_biotype/awk/extract_from_bed4_feature-E12/bedtools/makewindows_w-200_i-src/ChromHMM/LearnModel_test11_numstates-20_assembly-mm10/meta_20_segments/00_peak_anno_diagrams.pdf
#    """
#    input:
#        gtftk="opt/miniconda/envs/gtftk_dev/bin/gtftk",
#        chrominfo=input_chrominfo,
#        gtf=input_gtf,
#        bed="out/{filler}.bed"
#    output:
#        pdf="out/gtftk/peak_anno_chrominfo-{chrominfo_id}_gtf-{gtf_id}_more-keys-{more_keys}/{filler}/00_peak_anno_diagrams.pdf"
#    params:
#        outputdir="out/gtftk/peak_anno_chrominfo-{chrominfo_id}_gtf-{gtf_id}_more-keys-{more_keys}/{filler}"
#    threads:
#        1
#    wildcard_constraints:
#        gtf_id="[a-zA-Z0-9-]+",
#        chrominfo_id="[a-zA-Z0-9-]+",
#        more_keys="gene_biotype"
#    shell:
#        """
#        # Needed to get bedtools in the PATH:
#        export PATH="{WDIR}/opt/miniconda/envs/gtftk_dev/bin/:$PATH"
#
#        cd {params.outputdir}
#
#        {WDIR}/{input.gtftk} peak_anno \
#            --inputfile {WDIR}/{input.gtf} \
#            --more-keys {wildcards.more_keys} \
#            --outputdir . \
#            --chrom-info {WDIR}/{input.chrominfo} \
#            --peak-file {WDIR}/{input.bed} \
#            --no-date
#        """
#
#rule gtftk_peak_anno_chrominfo_gtf_more_keys_no_basic_feature:
#    """
#    Created: 
#        2017-06-26 15:47:51
#    Modified:
#    Aim:
#        Test --more-key feature.
#    Test:
#        out/gtftk/peak_anno_chrominfo-mm10-main-chr_gtf-GRCm38-main-chr_more-keys-gene_biotype_no-basic-feature/awk/extract_from_bed4_feature-E12/bedtools/makewindows_w-200_i-src/ChromHMM/LearnModel_test11_numstates-20_assembly-mm10/meta_20_segments/00_peak_anno_diagrams.pdf
#    """
#    input:
#        gtftk="opt/miniconda/envs/gtftk_dev/bin/gtftk",
#        chrominfo=input_chrominfo,
#        gtf=input_gtf,
#        bed="out/{filler}.bed"
#    output:
#        pdf="out/gtftk/peak_anno_chrominfo-{chrominfo_id}_gtf-{gtf_id}_more-keys-{more_keys}_no-basic-feature/{filler}/00_peak_anno_diagrams.pdf"
#    params:
#        outputdir="out/gtftk/peak_anno_chrominfo-{chrominfo_id}_gtf-{gtf_id}_more-keys-{more_keys}_no-basic-feature/{filler}"
#    threads:
#        1
#    wildcard_constraints:
#        gtf_id="[a-zA-Z0-9-]+",
#        chrominfo_id="[a-zA-Z0-9-]+",
#        more_keys="gene_biotype"
#    shell:
#        """
#        # Needed to get bedtools in the PATH:
#        export PATH="{WDIR}/opt/miniconda/envs/gtftk_dev/bin/:$PATH"
#
#        cd {params.outputdir}
#
#        {WDIR}/{input.gtftk} peak_anno \
#            --inputfile {WDIR}/{input.gtf} \
#            --more-keys {wildcards.more_keys} \
#            --outputdir . \
#            --chrom-info {WDIR}/{input.chrominfo} \
#            --no-basic-feature \
#            --peak-file {WDIR}/{input.bed} \
#            --no-date
#        """
#
#
