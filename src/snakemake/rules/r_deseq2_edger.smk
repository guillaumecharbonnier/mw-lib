"""
Important note:
    Have a look at dear.rules for differential expression analysis with DESeq2, EdgeR and Limma-voom.
"""

rule r_deseq2_edger:
    """
    Created:
        2019-02-05 16:58:37
    Aim:
        Deseq2 Edger
    Test:
        out/r/deseq2_edger/subread/featureCounts_-O_-t_exon_-g_gene_id_gtf-GRCh38-ensembl_bam-hg38-Casero2016-thy3-thy4/DESeq2_diagnostic_MA.png

    Check:
    git@github.com:wdecoster/DEA.R.git
        Script to automate differential expression analysis using DESeq2, edgeR or limma-voom


    """
    input:
        script = "src/r/deseq2_edger.R",
        #tsv="out/{filler}.tsv"
        tsv    = "out/{filler}.tsv"
    output:
        done   = "out/r/deseq2_edger/{filler}/done",
        ma     = "out/r/deseq2_edger/{filler}/DESeq2_diagnostic_MA.png",
        disp   = "out/r/deseq2_edger/{filler}/DESeq2_diagnostic_disp.png"
    params:
        outdir = "out/r/deseq2_edger/{filler}"
    conda:
        "../envs/r_deseq2_deger.R"
    shell:
        """
        {input.script} \
            --input_file {input.tsv} \
            --outdir {params.outdir} \
            --skip 2 \
            --class1 out/star/pe_GRCm38_outFilterMultimapNmax-1000/sickle/pe_-t_sanger_-q_20/merge_lanes/run176/RNA-R-H2AL2-WT-Rep1.bam,out/star/pe_GRCm38_outFilterMultimapNmax-1000/sickle/pe_-t_sanger_-q_20/merge_lanes/run176/RNA-R-H2AL2-WT-Rep2.bam,out/star/pe_GRCm38_outFilterMultimapNmax-1000/sickle/pe_-t_sanger_-q_20/merge_lanes/run176/RNA-R-H2AL2-WT-Rep3.bam \
            --class2 out/star/pe_GRCm38_outFilterMultimapNmax-1000/sickle/pe_-t_sanger_-q_20/merge_lanes/run176/RNA-R-H2AL2-KO-Rep1.bam,out/star/pe_GRCm38_outFilterMultimapNmax-1000/sickle/pe_-t_sanger_-q_20/merge_lanes/run176/RNA-R-H2AL2-KO-Rep2.bam,out/star/pe_GRCm38_outFilterMultimapNmax-1000/sickle/pe_-t_sanger_-q_20/merge_lanes/run176/RNA-R-H2AL2-KO-Rep3.bam

        touch {output.done}
        """




#rule r_deseq2_edger_featureCounts_gtfmerge-attr-GRCm38_texon_gmerge_gene_id_name_h2al2_rnaseq:
rule r_deseq2_edger_featureCounts_gtfmerge_attr_GRCm38_texon_gmerge_gene_id_name_h2al2_rnaseq:
    """
    Created:
        2017-02-14 15:49:38
    Modified:
        2017-05-02 10:44:44 - Generalization is removed because it is hard to write a simple rule which handle correctly classes.
    Aim:
        differential expression analysis using Denis' R script as basis.
    Mail Denis, important to consider for analysis:
        Hi,
        Lan Dao has performed a RNA-Seq analysis that I was in charge to analyse. I ran a DESeq2 analysis and was unable to find the knocked-down gene in the list of differentially expressed genes (while EdgeR called it significant). I googled a little bit to get some hints about this strange result. I end up with this (section 3.6 of DESeq2 manual*):

        RNA-seq data sometimes contain isolated instances of very large counts that are apparently unrelated to the experimental or study design, and which may be considered outliers... 
        Cook’s distance is intended to indicate an outlier count. The Cook’s distances are stored as a matrix available in assays(dds)[["cooks"]]....
        This filtering can be turned off with results(dds, cooksCutoff=FALSE).

        So I turned off the filter and ended up with my gene of interest whose adjusted p-value was ranked 1 after sorting (2.837314e-204). Then Lan and Salva were very happy (Hurrah! Hurrah!). Note that it also recover other genes that may be of interest.

        Questions:

           - Do you call that a bug ?
           - Are you a DESeq2 addict ?
           - Do you plan to use DESeq2 in the future ?
           - Would you recommend DESeq2 to your collaborators ?
           - Are statistics in DESeq2 and cooking related in some way ?
                        
            Ciao 

            *https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf

            The normalized counts for the gene of interest (XXX)

            K562_WT_rep1 14.35934553
            K562_WT_rep2 14.36950318
            K562_WT_rep3 14.69938411
            K562_XXX_clone105 3.054998913
            K562_XXX_clone109 2.610748962
                            
    Test:
        input:
            out/featureCounts/gtfmerge-attr-GRCm38_texon_gmerge_gene_id_name_h2al2_rnaseq.tsv
        output:
            out/r/deseq2_edger/featureCounts/gtfmerge-attr-GRCm38_texon_gmerge_gene_id_name_h2al2_rnaseq/done
    """
    input:
        rscript="opt/miniconda/envs/deseq2_edger/bin/Rscript",
        script="src/r/deseq2_edger.R",
        #tsv="out/{filler}.tsv"
        tsv="out/featureCounts/gtfmerge-attr-GRCm38_texon_gmerge_gene_id_name_h2al2_rnaseq.tsv"
    output:
        done="out/r/deseq2_edger/featureCounts/gtfmerge-attr-GRCm38_texon_gmerge_gene_id_name_h2al2_rnaseq/done",
        ma="out/r/deseq2_edger/featureCounts/gtfmerge-attr-GRCm38_texon_gmerge_gene_id_name_h2al2_rnaseq/DESeq2_diagnostic_MA.png",
        disp="out/r/deseq2_edger/featureCounts/gtfmerge-attr-GRCm38_texon_gmerge_gene_id_name_h2al2_rnaseq/DESeq2_diagnostic_disp.png"
    params:
        outdir="out/r/deseq2_edger/featureCounts/gtfmerge-attr-GRCm38_texon_gmerge_gene_id_name_h2al2_rnaseq"
    shell:
        """
        {input.rscript} {input.script} \
            --input_file {input.tsv} \
            --outdir {params.outdir} \
            --skip 2 \
            --class1 out/star/pe_GRCm38_outFilterMultimapNmax-1000/sickle/pe_-t_sanger_-q_20/merge_lanes/run176/RNA-R-H2AL2-WT-Rep1.bam,out/star/pe_GRCm38_outFilterMultimapNmax-1000/sickle/pe_-t_sanger_-q_20/merge_lanes/run176/RNA-R-H2AL2-WT-Rep2.bam,out/star/pe_GRCm38_outFilterMultimapNmax-1000/sickle/pe_-t_sanger_-q_20/merge_lanes/run176/RNA-R-H2AL2-WT-Rep3.bam \
            --class2 out/star/pe_GRCm38_outFilterMultimapNmax-1000/sickle/pe_-t_sanger_-q_20/merge_lanes/run176/RNA-R-H2AL2-KO-Rep1.bam,out/star/pe_GRCm38_outFilterMultimapNmax-1000/sickle/pe_-t_sanger_-q_20/merge_lanes/run176/RNA-R-H2AL2-KO-Rep2.bam,out/star/pe_GRCm38_outFilterMultimapNmax-1000/sickle/pe_-t_sanger_-q_20/merge_lanes/run176/RNA-R-H2AL2-KO-Rep3.bam

        touch {output.done}
        """

rule r_deseq2_edger_r_custom_table_post_blastn_counts_2017_04_24_wt_versus_ko:
    """
    Created:
        2017-02-14 15:49:38
    Modified:
        2017-05-02 10:44:44 - Generalization is removed because it is hard to write a simple rule which handle correctly classes.
    Aim:
        differential expression analysis using Denis R script as basis.
    Test:
        input:
            out/r/custom_table_post_blastn/counts_2017-04-24.tsv
        output:
            out/r/deseq2_edger/r/custom_table_post_blastn/counts_2017-04-24/WT-versus-KO-C-strandplus-pair1/done
            out/r/deseq2_edger/r/custom_table_post_blastn/counts_2017-04-24/WT-versus-KO-P-strandplus-pair1/done
            out/r/deseq2_edger/r/custom_table_post_blastn/counts_2017-04-24/WT-versus-KO-R-strandplus-pair1/done
            out/r/deseq2_edger/r/custom_table_post_blastn/counts_2017-04-24/WT-versus-KO-C-strandminus-pair1/done
            out/r/deseq2_edger/r/custom_table_post_blastn/counts_2017-04-24/WT-versus-KO-P-strandminus-pair1/done
            out/r/deseq2_edger/r/custom_table_post_blastn/counts_2017-04-24/WT-versus-KO-R-strandminus-pair1/done
    """
    input:
        rscript="opt/miniconda/envs/deseq2_edger/bin/Rscript",
        script="src/r/deseq2_edger.R",
        #tsv="out/{filler}.tsv"
        tsv="out/r/custom_table_post_blastn/counts_2017-04-24.tsv"
    output:
        ma_adj="out/r/deseq2_edger/r/custom_table_post_blastn/counts_2017-04-24/WT-versus-KO-{stage}-strand{strand}-pair{pair}/DESeq2_diagnostic_MA_adjusted.pdf",
        ma="out/r/deseq2_edger/r/custom_table_post_blastn/counts_2017-04-24/WT-versus-KO-{stage}-strand{strand}-pair{pair}/DESeq2_diagnostic_MA.pdf",
        disp="out/r/deseq2_edger/r/custom_table_post_blastn/counts_2017-04-24/WT-versus-KO-{stage}-strand{strand}-pair{pair}/DESeq2_diagnostic_disp.pdf"
    params:
        outdir="out/r/deseq2_edger/r/custom_table_post_blastn/counts_2017-04-24/WT-versus-KO-{stage}-strand{strand}-pair{pair}"
    wildcard_constraints:
        stage="P|R|C|S",
        strand="plus|minus",
        pair="1|2"
    shell:
        """
        if [ '{wildcards.stage}' = 'S' ]
        then
            CLASS1="RNA-{wildcards.stage}-H2AL2-WT-Rep1-strand{wildcards.strand}-pair{wildcards.pair},RNA-{wildcards.stage}-H2AL2-WT-Rep2-strand{wildcards.strand}-pair{wildcards.pair}"
            CLASS2="RNA-{wildcards.stage}-H2AL2-KO-Rep1-strand{wildcards.strand}-pair{wildcards.pair}"
        else
            CLASS1="RNA-{wildcards.stage}-H2AL2-WT-Rep1-strand{wildcards.strand}-pair{wildcards.pair},RNA-{wildcards.stage}-H2AL2-WT-Rep2-strand{wildcards.strand}-pair{wildcards.pair},RNA-{wildcards.stage}-H2AL2-WT-Rep3-strand{wildcards.strand}-pair{wildcards.pair}"
            CLASS2="RNA-{wildcards.stage}-H2AL2-KO-Rep1-strand{wildcards.strand}-pair{wildcards.pair},RNA-{wildcards.stage}-H2AL2-KO-Rep2-strand{wildcards.strand}-pair{wildcards.pair},RNA-{wildcards.stage}-H2AL2-KO-Rep3-strand{wildcards.strand}-pair{wildcards.pair}"
        fi
            

        {input.rscript} {input.script} \
            --input_file {input.tsv} \
            --outdir {params.outdir} \
            --skip 0 \
            --title H2AL2-WT-versus-KO-{wildcards.stage}-strand{wildcards.strand}-pair{wildcards.pair} \
            --class1 $CLASS1 \
            --class2 $CLASS2
        """

rule r_deseq2_edger_r_custom_table_post_blastn_counts_2017_04_24_stage_versus_stage:
    """
    Created:
        2017-02-14 15:49:38
    Modified:
        2017-05-02 10:44:44 - Generalization is removed because it is hard to write a simple rule which handle correctly classes.
    Aim:
        differential expression analysis using Denis R script as basis.
    """
    input:
        rscript="opt/miniconda/envs/deseq2_edger/bin/Rscript",
        script="src/r/deseq2_edger.R",
        #tsv="out/{filler}.tsv"
        tsv="out/r/custom_table_post_blastn/counts_2017-04-24.tsv"
    output:
        ma_adj="out/r/deseq2_edger/r/custom_table_post_blastn/counts_2017-04-24/stages-{stage1}-versus-{stage2}-strand{strand}-pair{pair}/DESeq2_diagnostic_MA_adjusted.pdf",
        ma="out/r/deseq2_edger/r/custom_table_post_blastn/counts_2017-04-24/stages-{stage1}-versus-{stage2}-strand{strand}-pair{pair}/DESeq2_diagnostic_MA.pdf",
        disp="out/r/deseq2_edger/r/custom_table_post_blastn/counts_2017-04-24/stages-{stage1}-versus-{stage2}-strand{strand}-pair{pair}/DESeq2_diagnostic_disp.pdf"
    params:
        outdir="out/r/deseq2_edger/r/custom_table_post_blastn/counts_2017-04-24/stages-{stage1}-versus-{stage2}-strand{strand}-pair{pair}"
    wildcard_constraints:
        stage1="P|R|C",
        stage2="P|R|C|S",
        strand="plus|minus",
        pair="1|2"
    shell:
        """
        if [ '{wildcards.stage2}' = 'S' ]
        then
            CLASS2="RNA-{wildcards.stage2}-H2AL2-WT-Rep1-strand{wildcards.strand}-pair{wildcards.pair},RNA-{wildcards.stage2}-H2AL2-WT-Rep2-strand{wildcards.strand}-pair{wildcards.pair},RNA-{wildcards.stage2}-H2AL2-KO-Rep1-strand{wildcards.strand}-pair{wildcards.pair}"
        else
            CLASS2="RNA-{wildcards.stage2}-H2AL2-WT-Rep1-strand{wildcards.strand}-pair{wildcards.pair},RNA-{wildcards.stage2}-H2AL2-WT-Rep2-strand{wildcards.strand}-pair{wildcards.pair},RNA-{wildcards.stage2}-H2AL2-WT-Rep3-strand{wildcards.strand}-pair{wildcards.pair},RNA-{wildcards.stage2}-H2AL2-KO-Rep1-strand{wildcards.strand}-pair{wildcards.pair},RNA-{wildcards.stage2}-H2AL2-KO-Rep2-strand{wildcards.strand}-pair{wildcards.pair},RNA-{wildcards.stage2}-H2AL2-KO-Rep3-strand{wildcards.strand}-pair{wildcards.pair}"
        fi
        
        CLASS1="RNA-{wildcards.stage1}-H2AL2-WT-Rep1-strand{wildcards.strand}-pair{wildcards.pair},RNA-{wildcards.stage1}-H2AL2-WT-Rep2-strand{wildcards.strand}-pair{wildcards.pair},RNA-{wildcards.stage1}-H2AL2-WT-Rep3-strand{wildcards.strand}-pair{wildcards.pair},RNA-{wildcards.stage1}-H2AL2-KO-Rep1-strand{wildcards.strand}-pair{wildcards.pair},RNA-{wildcards.stage1}-H2AL2-KO-Rep2-strand{wildcards.strand}-pair{wildcards.pair},RNA-{wildcards.stage1}-H2AL2-KO-Rep3-strand{wildcards.strand}-pair{wildcards.pair}"
        
        {input.rscript} {input.script} \
            --input_file {input.tsv} \
            --outdir {params.outdir} \
            --skip 0 \
            --title H2AL2-stages-{wildcards.stage1}-versus-{wildcards.stage2}-strand{wildcards.strand}-pair{wildcards.pair} \
            --class1 $CLASS1 \
            --class2 $CLASS2
        """

rule tmp_deseq2:
    input:
        "out/r/deseq2_edger/r/custom_table_post_blastn/counts_2017-04-24/WT-versus-KO-S-strandplus-pair1/DESeq2_diagnostic_MA_adjusted.pdf",
        "out/r/deseq2_edger/r/custom_table_post_blastn/counts_2017-04-24/WT-versus-KO-C-strandplus-pair1/DESeq2_diagnostic_MA_adjusted.pdf",
        "out/r/deseq2_edger/r/custom_table_post_blastn/counts_2017-04-24/WT-versus-KO-P-strandplus-pair1/DESeq2_diagnostic_MA_adjusted.pdf",
        "out/r/deseq2_edger/r/custom_table_post_blastn/counts_2017-04-24/WT-versus-KO-R-strandplus-pair1/DESeq2_diagnostic_MA_adjusted.pdf",
        "out/r/deseq2_edger/r/custom_table_post_blastn/counts_2017-04-24/WT-versus-KO-C-strandminus-pair1/DESeq2_diagnostic_MA_adjusted.pdf",
        "out/r/deseq2_edger/r/custom_table_post_blastn/counts_2017-04-24/WT-versus-KO-P-strandminus-pair1/DESeq2_diagnostic_MA_adjusted.pdf",
        "out/r/deseq2_edger/r/custom_table_post_blastn/counts_2017-04-24/WT-versus-KO-S-strandminus-pair1/DESeq2_diagnostic_MA_adjusted.pdf",
        "out/r/deseq2_edger/r/custom_table_post_blastn/counts_2017-04-24/WT-versus-KO-R-strandminus-pair1/DESeq2_diagnostic_MA_adjusted.pdf",
        expand("out/r/deseq2_edger/r/custom_table_post_blastn/counts_2017-04-24/stages-{stage1}-versus-{stage2}-strand{strand}-pair{pair}/DESeq2_diagnostic_MA_adjusted.pdf",
            stage1=["P","R","C"],
            stage2=["P","R","C","S"],
            strand=["plus", "minus"],
            pair=["1"])

