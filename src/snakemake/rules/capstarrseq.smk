"""
Created:
    2017-04-11 11:21:34
Aim:
    Workflow for capstarseq analysis
    Merge my workflow for Grenoble and the Capstarrseq for Salva to ease development of rules.
Taken from:
    salva/fq_to_bam/code/snakemake/workflows/capstarseq_2017_01_05.snakefile
Notes:
    Capstarrseq analysis for run170 and rerun of run150 without FPKM filtering.
    Extension: 314
    Mouse mm9
    No filter on input FPKM
    Previous CapStarrseq_mSilencers (run 150)
    /gpfs/tgml/reads/fastq/Run_150_NS500-057_12.09.2016_SS_BSoS
    
    DHS_PGK_DNA_library (change name to= mT_DHS_PGK_P5424_rep1)
    Input input_plasmid_library (change name to= mT_DHS_PGK_input)
    
    Run 170:
    MT-DHS-PGK-p5424_rep2
    Input: input_plasmid_library (mT_DHS_PGK_input) from run 150
    
    MT-DHS-EFIA-p5424_rep1
    MT-DHS-EFIA-p5424_rep2
    MT-DHS-EFIA-p5424_rep3
    Input: mT-DHS-EFIA-input
    

"""

EXT_CAPSTARSEQ = "314"
EXT_ATAC = "180"
EXT_CHIPSEQ = "300"
EXT_IGMM = "437"



ID_BAM_TO_BED = "|".join([
    "bedtools/bamtobed|awk/extend_reads_[0-9]+/bedtools/bamtobed/samtools/sort_-n",
    "awk/convert_bedpe_to_bed6_insert_size/bedtools/bamtobed_bedpe/samtools/sort_-n",
    "awk/extend_reads_[0-9]+/awk/keep_first_mate_for_pe_bedtools_bamtobed/bedtools/bamtobed/sort_-n",
    "awk/convert_bedpe_to_bed6_insert_size/bedtools/bamtobed_-bedpe/samtools/sort_-n"])

CRM_TYPE = "mTDHS|hProm|hProm_posEprom|IGMM|hSE_RPMI_JURKAT_DND41_hg19_Alex"

#rule capstarrseq_coverage_fpkm_input_v2:
#    input:
#        
#

rule capstarrseq_coverage_fpkm_input:
    """
    Modified:
        2018-01-30 13:34:24 - FPKM calculation does not rely on flagstat number of mapped reads anymore. It seems to be a better solution with just a wc -l on input bed for paired-end data merged with single-end.    
    Test:
        out/15_mTDHS/awk/extend_reads_314/bedtools/bamtobed/samtools/sam_to_bam/bowtie2/pe_mm9/sickle/pe_-t_sanger_-q_20/gunzip/ln/rename_run178_tgml/BAC_mT_DHS_PGK_p5424_rep1.filtered_CRMs.bed
    Note:
        I don't understand why flagstat is used here. Maybe juste a wc -l of the input bed would be more appropriate. Especially because flagstat is not appropriate when only mate 1 of paired-end reeds are kepts:
        Number of reads based on flagstat:
        47449907
        39311821 out/sort/coordinates_bed/awk/extend_reads_314/awk/keep_first_mate_for_pe_bedtools_bamtobed/bedtools/bamtobed/ln/sst_exp/CapStarr_155_170_178_204/mm9/mT_DHS_PGK_input.bed

    """
    input:
        #bed_reads="out/{id_bam_to_bed}/{id}.bed",
        #bed_crms="inp/bed/crms/{crm_type}.bed",
        bed_reads="out/sort/coordinates_bed/{id_bam_to_bed}/{id}.bed",
        bed_crms="out/sort/coordinates_bed/ln/updir/mw/inp/bed/crms/{crm_type}.bed",
        flagstat="out/{id}.flagstat.txt"
    output:
        tsv_cov_unfiltered="out/capstarrseq/coverage_fpkm_input_{crm_type}/{id_bam_to_bed}/{id}.coverage_unfiltered.tsv",
        tsv_fpkm="out/capstarrseq/coverage_fpkm_input_{crm_type}/{id_bam_to_bed}/{id}.FPKM.tsv",
        tsv_fpkm_unfiltered="out/capstarrseq/coverage_fpkm_input_{crm_type}/{id_bam_to_bed}/{id}.FPKM_unfiltered.tsv",
        bed_crms="out/capstarrseq/coverage_fpkm_input_{crm_type}/{id_bam_to_bed}/{id}.filtered_CRMs.bed"
    params:
        fpkm_threshold='1'
    wildcard_constraints:
        id_bam_to_bed = ID_BAM_TO_BED, 
        crm_type = CRM_TYPE
    conda:
        "../envs/capstarrseq.yaml"
    shell:
        """
        # 1. Coverage
        # Old way consume too much RAM.
        #bedtools coverage -a {input.bed_crms} -b {input.bed_reads} > {output.tsv_cov_unfiltered}
        # Strange because it still consumes a lot of memory with 'sorted'.
        bedtools coverage -sorted -a {input.bed_crms} -b {input.bed_reads} > {output.tsv_cov_unfiltered}
        #cp {output.tsv_cov_unfiltered} tmp/debug_tsv_cov_unfiltered.tsv
        
        # 2. Coverage to FPKM
        nb_mReads_flagstat=`grep "mapped (" {input.flagstat} | awk '{{print $1}}'`
        nb_mReads=`wc -l {input.bed_reads} | cut -f1 -d ' '`
        echo "Number of reads based on flagstat: $nb_mReads_flagstat"
        echo "Number of reads based on wc -l   : $nb_mReads"
        echo "debug: check the difference between the two values for number of reads"

        awk -v NB_READS=$nb_mReads '{{ fpkm = ($(NF-3) / (($(NF-1))/1000)) / (NB_READS/1000000) ; print $1"\\t"$2"\\t"$3"\\t"$4"\\t"fpkm }}' {output.tsv_cov_unfiltered} > {output.tsv_fpkm_unfiltered}
        #cp {output.tsv_fpkm_unfiltered} tmp/debug_capstarrseq_bedtools_coverage_fpkm_input.tsv
        
        Rscript -e '
            print("3. CRM filtering based on FPKM threshold (input FPKM < 1)")
            fpkm_all <- read.table(
                "{output.tsv_fpkm_unfiltered}",
                stringsAsFactors = FALSE
            )
            idx <- which(fpkm_all[,5] >= {params.fpkm_threshold})
            fpkm_all <- fpkm_all[idx,]
            write.table(
                fpkm_all,
                file = "{output.tsv_fpkm}",
                quote = FALSE,
                row.names = FALSE,
                col.names = FALSE,
                sep = "\\t"
            )

            print("4. Producing filtered CRM file")
            crms_all <- read.table(
                "{input.bed_crms}",
                stringsAsFactors = FALSE,
                sep = "\\t"
            )
            #fpkm_all <- read.table("{output.tsv_fpkm}", stringsAsFactors=F)
            print("read.table done")
            str_crms_all <- apply(
                crms_all[,1:3],
                1,
                paste,
                collapse="_"
            )
            str_fpkm_all <- apply(
                fpkm_all[,1:3],
                1,
                paste,
                collapse = "_"
            )
            idx <- match(
                str_fpkm_all,
                str_crms_all
            )
            write.table(
                crms_all[idx,],
                file = "{output.bed_crms}",
                quote = FALSE,
                row.names = FALSE,
                col.names = FALSE,
                sep = "\\t"
            )
        '
        """

rule capstarrseq_coverage_fpkm_sample:
    """
    Modified:
        2018-01-30 13:34:24 - FPKM calculation does not rely on flagstat number of mapped reads anymore. It seems to be a better solution with just a wc -l on input bed for paired-end data merged with single-end.
    """
    input:
        bed_reads="out/sort/coordinates_bed/{id_bam_to_bed}/{id}/{id_sample}.bed",
        bed_crms="out/capstarrseq/coverage_fpkm_input_{crm_type}/{id_bam_to_bed}/{id}/{id_input}.filtered_CRMs.bed",
        flagstat="out/{id}/{id_sample}.flagstat.txt"
    output:
        tsv_cov  = "out/capstarrseq/coverage_fpkm_sample_{crm_type}/{id_bam_to_bed}/{id}/{id_sample}_on_{id_input}_filtered_CRMs.coverage.tsv",
        tsv_fpkm = "out/capstarrseq/coverage_fpkm_sample_{crm_type}/{id_bam_to_bed}/{id}/{id_sample}_on_{id_input}_filtered_CRMs.FPKM.tsv"
    wildcard_constraints:
        id_bam_to_bed = ID_BAM_TO_BED,
        crm_type = CRM_TYPE
    conda:
        "../envs/capstarrseq.yaml"
    shell:
        """
        # 1. Coverage
        bedtools coverage -sorted -a {input.bed_crms} -b {input.bed_reads} > {output.tsv_cov}

        # 2. Coverage to FPKM
        nb_mReads_flagstat=`grep "mapped (" {input.flagstat} | awk '{{print $1}}'`
        nb_mReads=`wc -l {input.bed_reads} | cut -f1 -d ' '`
        echo "Number of reads based on flagstat: $nb_mReads_flagstat"
        echo "Number of reads based on wc -l   : $nb_mReads"
        echo "debug: check the difference between the two values for number of reads"

        awk -v NB_READS=$nb_mReads '{{ fpkm = ($(NF-3) / (($(NF-1))/1000)) / (NB_READS/1000000) ; print $1"\\t"$2"\\t"$3"\\t"$4"\\t"fpkm }}' {output.tsv_cov} > {output.tsv_fpkm}
        """

rule capstarrseq_fold_change:
    """
    Created:
        2016 - From Aurelien's workflow.
    Modified:
        2018-01-31 15:46:20 - Added +1 count to both fpkm in order to prevent infinite value issue.
    Aim:
        Computation of fold change in samples.
    """
    input:
        tsv_fpkm_sample = "out/capstarrseq/coverage_fpkm_sample_{crm_type}/{id}/{id_sample}_on_{id_input}_filtered_CRMs.FPKM.tsv",
        tsv_fpkm_input  = "out/capstarrseq/coverage_fpkm_input_{crm_type}/{id}/{id_input}.FPKM.tsv"
    output:
        fc = "out/capstarrseq/fold_change_{crm_type}/{id}/{id_sample}_over_{id_input}.foldChange.tsv"
    wildcard_constraints:
        #id_bam_to_bed="bedtools/bamtobed|awk/extend_reads_[0-9]+/bedtools/bamtobed",
        crm_type = CRM_TYPE
    conda:
        "../envs/capstarrseq.yaml"
    shell:
        """
        Rscript -e '
            fpkm_input <- read.table(
                "{input.tsv_fpkm_input}",
                stringsAsFactors = FALSE
            )
            fpkm_sample <- read.table(
                "{input.tsv_fpkm_sample}",
                stringsAsFactors = FALSE
            )
            fc <- (fpkm_sample[,5] + 1) / (fpkm_input[,5] + 1)
            dat <- data.frame(fpkm_sample[,1:4], fc)
            write.table(
                dat,
                file = "{output.fc}",
                quote = FALSE,
                row.names = FALSE,
                col.names = FALSE,
                sep = "\\t"
            )
        '
        """

rule capstarrseq_grouping_crms:
    """
    Created:
        2016 - From Aurelien's workflow.
    Modified:
        2018-01-29 16:39:47
    Aim:
        Creation des groupes de CRMs (inactifs, actifs) sur la base des Fold Change et des categories
    Test:
        out/capstarrseq/grouping_crms/capstarrseq/fold_change_mTDHS/awk/extend_reads_314/bedtools/bamtobed/ln/sst_exp/CapStarr_155_170_178_204/mm9/mT_DHS_PGK_P5424_stimulated_rep3_over_mT_DHS_PGK_input.inflexionPointGroups.tsv
    """
    input:
        fc = "out/{id}/{id_sample}_over_{id_input}.foldChange.tsv"
    output:
        groups = "out/capstarrseq/grouping_crms/{id}/{id_sample}_over_{id_input}.inflexionPointGroups.tsv",
        pdf = "out/capstarrseq/grouping_crms/{id}/{id_sample}_over_{id_input}.inflexionPointGroups.pdf"
    params:
        outdir = "out/capstarrseq/grouping_crms/{id}/",
        th_negative_1 = "0.05",
        th_negative_2 = "0.01"
    conda:
        "../envs/capstarrseq.yaml"
    shell:
        """
        Rscript -e '
            dat <- read.table(
                "{input.fc}",
                stringsAsFactors = FALSE
            )
        
            ## CALCUL DU SEUIL EN FONCTION DU FDR
            for (fdr in c({params.th_negative_1},{params.th_negative_2})) {{
                # Determination du seuil
                idx <- which(dat[,4] == "Negative")
                # Skipping FDR computation if no Negative set is provided.
                if (length(idx)!=0) {{
                    P <- ecdf(dat[idx,5])
                    th_FoldChange <- quantile(P,probs=1-fdr)

                    # identification des regions actives/inactives
                    idx <- which(dat[,5] >= th_FoldChange)
                    groups <- rep(
                        "Inactive",
                        nrow(dat)
                    )
                    groups[idx] <- "Active"
                    groups[which(is.na(dat[,5]))] <- "NA"
                    
                    pdf(
                        paste(
                            "{params.outdir}/Groups.FDR=",
                            fdr,
                            ".pdf",
                            sep=""
                        )
                    )
                    par(mfrow=c(2,2))
                    
                    #- boxplot en fonction des categories
                    categories <- unique(dat[,4])
                    fc <- list() ; lfc <- list()
                    for (category in categories) {{
                        idx <- which(dat[,4] == category)
                        fc[[category]] <- dat[idx,5]
                        lfc[[category]] <- log2(dat[idx,5])
                    }}
                    colors <- rep("lightgrey",length(categories))
                    if ("Random" %in% categories) {{ colors[which(categories == "Random")] <- "darkblue" }}
                    if ("PosEpromoter" %in% categories) {{ colors[which(categories == "PosEpromoter")] <- "darkgreen" }}
                    if ("Positive" %in% categories) {{ colors[which(categories == "Positive")] <- "darkgreen" }}
                    if ("Negative" %in% categories) {{ colors[which(categories == "Negative")] <- "darkred" }}
                    boxplot(
                        fc,
                        pch = 20,
                        col=colors,
                        main = "Fold Change of categories",
                        ylab = "Fold Change",
                        las = 2
                    )
                    text(
                        length(categories),
                        0.9*max(dat[,5]),
                        labels = sprintf(
                            "Threshold :\n%3.2f",
                            th_FoldChange
                        ),
                        col = "red"
                    )
                    boxplot(
                        lfc,
                        pch = 20,
                        col = colors,
                        main = "Fold Change of categories",
                        ylab = "Fold Change [log2]",
                        las = 2
                    )
                    abline(
                        h = log2(th_FoldChange),
                        col = "red",
                        lty = "dashed"
                    )

                    #- ranked genomic regions based on their Fold Change
                    plot(
                        sort(
                            dat[
                                -which(dat[,4]=="Random" | dat[,4]=="Negative"),
                                5
                            ]
                        ),
                        pch = 20,
                        main = "Activity of genomic regions\n(Random/Negative not included)",
                        ylab = "Fold Change",
                        xlab = "Ranked genomic regions"
                    )
                    idx <- which(dat[,5] >= th_FoldChange & dat[,4] != "Random" & dat[,4] != "Negative")
                    abline(
                        v = nrow(dat)-length(idx),
                        col = "red",
                        lty = "dashed"
                    )
                    text(
                        nrow(dat)-length(idx),
                        0.9*max(dat[,5]),
                        labels = length(idx),
                        pos = 2,
                        offset = 0,
                        col = "red"
                    )
                    idx <- which(
                        dat[,5] < th_FoldChange &
                        dat[,4] != "Random" &
                        dat[,4] != "Negative"
                    )
                    text(
                        length(idx)/2,
                        0.9*max(dat[,5]),
                        labels = length(idx),
                        col = "black"
                    )
            
                    #- info sur l"analyse
                    plot(
                        c(0, 1),
                        c(0, 1),
                        ann = FALSE,
                        bty = "n",
                        type = "n",
                        xaxt = "n",
                        yaxt = "n"
                    )
                    text(
                        x = 0.5,
                        y = 1,
                        "INFORMATION",
                        cex = 1.5,
                        pos = 1,
                        offset = 0,
                        col = "black"
                    )
                    text(
                        x = 0,
                        y = 0.70,
                        "{wildcards.id_sample}",
                        pos = 4,
                        offset = 0,
                        col = "black"
                    )
                    text(
                        x = 0,
                        y = 0.60,
                        paste(
                            "Method: FDR=",
                            fdr,
                            sep = ""
                        ),
                        pos = 4,
                        offset = 0,
                        col = "black"
                    )
                    dev.off()
                    
                    dat_groups <- data.frame(
                        dat[,1:4],
                        groups
                    )
                    write.table(
                        dat_groups,
                        file = paste(
                            "{params.outdir}/Groups.FDR=",
                            fdr,
                            ".grp",
                            sep = ""
                        ),
                        quote = FALSE,
                        row.names = FALSE,
                        col.names = FALSE,
                        sep = "\\t"
                    )
                }}
            }}
            
            
            ## Compute threshold base on inflexion point
            
            #---  code from ROSE tool to dertermine super-enhancers  ---#
            #--
            #This function calculates the cutoff by sliding a diagonal line and finding where it is tangential (or as close as possible)
            calculate_cutoff <- function(inputVector, drawPlot=TRUE,...){{
                inputVector <- sort(inputVector)
                inputVector[inputVector<0]<-0 #set those regions with more control than ranking equal to zero
                slope <- (max(inputVector)-min(inputVector))/length(inputVector) #This is the slope of the line we want to slide. This is the diagonal.
                xPt <- floor(optimize(numPts_below_line,lower=1,upper=length(inputVector),myVector= inputVector,slope=slope)$minimum) #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
                y_cutoff <- inputVector[xPt] #The y-value at this x point. This is our cutoff.
        
                if(drawPlot){{  #if TRUE, draw the plot
                    plot(1:length(inputVector), inputVector,type="l",...)
                    b <- y_cutoff-(slope* xPt)
                    abline(
                        v = xPt,
                        h = y_cutoff,
                        lty = 2,
                        col = 8
                    )
                    points(
                        xPt,
                        y_cutoff,
                        pch = 16,
                        cex = 0.9,
                        col = 2
                    )
                    abline(
                        coef = c(b,slope),
                        col = 2
                    )
                    title(
                        paste(
                            "x=",
                            xPt,
                            "\ny=",
                            signif(
                                y_cutoff,
                                3
                            ),
                            "\nFold over Median=",
                            signif(
                                y_cutoff/median(inputVector),
                                3
                            ),
                            "x\nFold over Mean=",
                            signif(
                                y_cutoff/mean(inputVector),
                                3
                            ),
                            "x",
                            sep=""
                        )
                    )
                    axis(
                        1,
                        sum(inputVector==0),
                        sum(inputVector==0),
                        col.axis = "pink",
                        col = "pink"
                    ) #Number of regions with zero signal
                }}
                return(
                    list(
                        absolute = y_cutoff,
                        overMedian = y_cutoff/median(inputVector),
                        overMean = y_cutoff/mean(inputVector)
                    )
                )
            }}

            #this is an accessory function, that determines the number of points below a diagnoal passing through [x,yPt]
            numPts_below_line <- function(myVector,slope,x){{
                yPt <- myVector[x]
                b <- yPt-(slope*x)
                xPts <- 1:length(myVector)
                return(sum(myVector<=(xPts*slope+b)))
            }}
            #--
            #---  code from ROSE tool to determine super-enhancers  ---#

            # determination du seuil
            idx <- which(
                dat[,4] != "Negative" &
                dat[,4] != "Random"
            )
            th_FoldChange <- calculate_cutoff(dat[idx,5], drawPlot=F)$absolute

            # identification des regions actives/inactives
            idx <- which(dat[,5] >= th_FoldChange)
            groups <- rep(
                "Inactive",
                nrow(dat)
            )
            groups[idx] <- "Active"
            groups[which(is.na(dat[,5]))] <- "NA"

            pdf("{output.pdf}")
            par(mfrow=c(2,2))

            #- boxplot en fonction des categories
            categories <- unique(dat[,4])
            fc <- list() ; lfc <- list()
            for (category in categories) {{
                idx <- which(dat[,4] == category)
                fc[[category]] <- dat[idx,5]
                lfc[[category]] <- log2(dat[idx,5])
            }}
            colors <- rep("lightgrey",length(categories))
            if ("Random" %in% categories) {{ colors[which(categories == "Random")] <- "darkblue" }}
            if ("PosEpromoter" %in% categories) {{ colors[which(categories == "PosEpromoter")] <- "darkgreen" }}
            if ("Positive" %in% categories) {{ colors[which(categories == "Positive")] <- "darkgreen" }}
            if ("Negative" %in% categories) {{ colors[which(categories == "Negative")] <- "darkred" }}
            boxplot(fc, pch=20, col=colors, main="Fold Change of categories", ylab="Fold Change", las=2)
            text(length(categories), 0.9*max(dat[,5]), labels=sprintf("Threshold :\n%3.2f",th_FoldChange), col="red")
            boxplot(lfc, pch=20, col=colors, main="Fold Change of categories", ylab="Fold Change [log2]", las=2)
            abline(h=log2(th_FoldChange), col="red", lty="dashed")
            
            #- ranked genomic regions based on their Fold Change
            plot(sort(dat[which(dat[,4]!="Random" & dat[,4]!="Negative"),5]), pch=20, main="Activity of genomic regions\n(Random/Negative not included)", ylab="Fold Change", xlab="Ranked genomic regions")

            idx <- which(dat[,5] >= th_FoldChange & dat[,4] != "Random" & dat[,4] != "Negative")
            abline(v=nrow(dat)-length(idx), col="red", lty="dashed")
            text(nrow(dat)-length(idx)/2, 0.9*max(dat[,5]), labels=length(idx), pos=2, offset=0, col="red")
            idx <- which(dat[,5] < th_FoldChange & dat[,4] != "Random" & dat[,4] != "Negative")
            text(length(idx)/2, 0.9*max(dat[,5]), labels=length(idx), col="black")
            
            print("debug_test4")
            
            #- info sur l"analyse
            plot(
                c(0, 1),
                c(0, 1),
                ann = FALSE,
                bty = "n",
                type = "n",
                xaxt = "n",
                yaxt = "n"
            )
            text(
                x = 0.5,
                y = 1,
                "INFORMATION",
                cex = 1.5,
                pos = 1,
                offset = 0,
                col = "black"
            )
            text(
                x = 0,
                y = 0.70,
                "{wildcards.id_sample}",
                pos = 4,
                offset = 0,
                col = "black"
            )
            text(
                x = 0,
                y = 0.60,
                "Method: Inflexion point",
                pos = 4,
                offset = 0,
                col = "black"
            )
            dev.off()

            dat <- data.frame(
                dat[,1:4],
                groups
            )
            write.table(
                dat,
                file = "{output.groups}",
                quote = FALSE,
                row.names = FALSE,
                col.names = FALSE,
                sep = "\\t"
            )
        '
        """

rule capstarrseq_merge_all_data:
    """
    Test:
        out/capstarrseq/coverage_fpkm_input_hProm_posEprom/awk/extend_reads_314/bedtools/bamtobed/samtools/sam_to_bam/bowtie2/se_hg19/sickle/se_-t_sanger_-q_20/gunzip/merge_lanes_nextseq500_single_end/ln/rename_run107_tgml/Input.filtered_CRMs.bed
        out/capstarrseq/merge_all_data_mTDHS/awk/extend_reads_314/bedtools/bamtobed/ln/sst_exp/CapStarr_155_170_178_204/mm9/mT_DHS_PGK_P5424_stimulated_rep3_over_mT_DHS_PGK_input.tsv
        out/capstarrseq/grouping_crms/capstarrseq/fold_change_mTDHS/awk/extend_reads_314/bedtools/bamtobed/ln/sst_exp/CapStarr_155_170_178_204/mm9/mT_DHS_PGK_P5424_stimulated_rep3_over_mT_DHS_PGK_input.inflexionPointGroups.tsv
    Nori 2019-09-05 13:54:50:
        out/capstarrseq/merge_all_data_hProm/awk/convert_bedpe_to_bed6_insert_size/bedtools/bamtobed_-bedpe/samtools/sort_-n/ln/alias/sst/all_samples/hg19/bam/CapStarr_K562_IFN_rep1_Seq2_over_Capstarr_K562_IFN_control.allData.tsv out/capstarrseq/merge_all_data_hProm/awk/convert_bedpe_to_bed6_insert_size/bedtools/bamtobed_-bedpe/samtools/sort_-n/ln/alias/sst/all_samples/hg19/bam/CapStarr_K562_IFN_rep2_Seq2_over_Capstarr_K562_IFN_control.allData.tsv

        out/capstarrseq/merge_all_data_hProm/awk/extend_reads_314/bedtools/bamtobed/samtools/merge/bam-hg19-CapStarr-K562-IFN-rep1-merged_over_RAJOUTERIDCONTROLMERGED.allData.tsv

    Ahmad 2019-10-15 00:34:22
        out/capstarrseq/merge_all_data_IGMM/awk/convert_bedpe_to_bed6_insert_size/bedtools/bamtobed_-bedpe/samtools/sort_-n/ln/alias/sst/all_samples/mm9/bam/capSTARR-seq_P5424_rep1_over_capSTARR-seq_P5424_Input.allData.tsv out/capstarrseq/merge_all_data_IGMM/awk/convert_bedpe_to_bed6_insert_size/bedtools/bamtobed_-bedpe/samtools/sort_-n/ln/alias/sst/all_samples/mm9/bam/capSTARR-seq_P5424_rep2_over_capSTARR-seq_P5424_Input.allData.tsv
    """
    input:
        bed_crm="out/capstarrseq/coverage_fpkm_input_{crm_type}/{id}/{id_input}.filtered_CRMs.bed",
        fpkm = "out/capstarrseq/coverage_fpkm_sample_{crm_type}/{id}/{id_sample}_on_{id_input}_filtered_CRMs.FPKM.tsv",
        fpkm_input = "out/capstarrseq/coverage_fpkm_input_{crm_type}/{id}/{id_input}.FPKM.tsv",
        fc = "out/capstarrseq/fold_change_{crm_type}/{id}/{id_sample}_over_{id_input}.foldChange.tsv",
        group = "out/capstarrseq/grouping_crms/capstarrseq/fold_change_{crm_type}/{id}/{id_sample}_over_{id_input}.inflexionPointGroups.tsv",
        pdf = "out/capstarrseq/grouping_crms/capstarrseq/fold_change_{crm_type}/{id}/{id_sample}_over_{id_input}.inflexionPointGroups.pdf",
    output:
        tsv = "out/capstarrseq/merge_all_data_{crm_type}/{id}/{id_sample}_over_{id_input}.allData.tsv",
        pdf = "out/capstarrseq/merge_all_data_{crm_type}/{id}/{id_sample}_over_{id_input}.inflexionPointGroups.pdf"
    conda:
        "../envs/capstarrseq.yaml"
    shell:
        """
        cp {input.pdf} {output.pdf}
        Rscript -e '
            bed_crm <- read.table(
                "{input.bed_crm}",
                stringsAsFactors = FALSE,
                sep = "\\t"
            )
            fpkm <- read.table(
                "{input.fpkm}",
                stringsAsFactors = FALSE
            )
            fpkm_input <- read.table(
                "{input.fpkm_input}",
                stringsAsFactors = FALSE
            )
            fc <- read.table(
                "{input.fc}",
                stringsAsFactors = FALSE
            )
            group <- read.table(
                "{input.group}",
                stringsAsFactors = FALSE
            )

            # Because some of my input bed are bed4 and some are bed6.
            if (dim(bed_crm)[2] == 4){{
                bed_crm[,5] <- "."
                bed_crm[,6] <- "NA"
                }}

            dat <- data.frame(
                fpkm,
                fpkm_input[,5],
                fc[,5],
                group[,5],
                bed_crm[,6]
            )
            #group_labels <- sapply(
            #    strsplit(
            #        group_files,
            #        "\\\."
            #    ),
            #    function(x){{
            #        nb = length(x) ;
            #        return(
            #            paste(
            #                "group_",
            #                x[nb-2],
            #                ".",
            #                x[nb-1],
            #                sep = ""
            #            )
            #        )
            #    }}
            #)
            colnames(dat) <- c(
                "chr",
                "start",
                "end",
                "category",
                "fpkm",
                "fpkm_input",
                "fold_change",
                "group_inflexionPoint",
                "gene_name"
            )
            
            #if (file.exists("data/CRMs/CRMs_geneAssociation.bed_notOrdered") == TRUE) {{
            #    genes <- read.table("data/CRMs/CRMs_geneAssociation.bed_notOrdered", stringsAsFactors=F)
            #    str_fpkm <- apply(fpkm[,1:3],1,paste,collapse="_")
            #    str_genes <- apply(genes[,1:3],1,paste,collapse="_")
            #    idx <- match(str_fpkm, str_genes)
            #    genes <- genes[idx,]
            #    dat <- data.frame(dat, genes=genes[,4], location=genes[,5])
            #}}

            write.table(
                dat,
                file = "{output.tsv}",
                quote = FALSE,
                row.names = FALSE,
                col.names = TRUE,
                sep="\\t"
            )
        '
        """
