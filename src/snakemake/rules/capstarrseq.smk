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
        # 1. Coverage
        # Old way consume too much RAM.
        #{input.bedtools} coverage -a {input.bed_crms} -b {input.bed_reads} > {output.tsv_cov_unfiltered}
        # Strange because it still consumes a lot of memory with 'sorted'.
        {input.bedtools} coverage -sorted -a {input.bed_crms} -b {input.bed_reads} > {output.tsv_cov_unfiltered}
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
    run:
        R("""
        fpkm_input <- read.table('{input.tsv_fpkm_input}', stringsAsFactors=F)
        fpkm_sample <- read.table('{input.tsv_fpkm_sample}', stringsAsFactors=F)
        fc <- (fpkm_sample[,5] + 1) / (fpkm_input[,5] + 1)
        dat <- data.frame(fpkm_sample[,1:4], fc)
        write.table(dat, file='{output.fc}', quote=F, row.names=F, col.names=F, sep='\t')
        """)

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
    run:
        R("""
        dat <- read.table('{input.fc}', stringsAsFactors=F)
    
        ## CALCUL DU SEUIL EN FONCTION D'UN FDR
        for (fdr in c({params.th_negative_1},{params.th_negative_2})) {{
            # Determination du seuil
            idx <- which(dat[,4] == 'Negative')
            # Skipping FDR computation if no Negative set is provided.
            if (length(idx)!=0) {{
                P <- ecdf(dat[idx,5])
                th_FoldChange <- quantile(P,probs=1-fdr)

                # identification des regions actives/inactives
                idx <- which(dat[,5] >= th_FoldChange)
                groups <- rep('Inactive',nrow(dat))
                groups[idx] <- 'Active'
                groups[which(is.na(dat[,5]))] <- 'NA'
                
                pdf(paste('{params.outdir}/Groups.FDR=',fdr,'.pdf',sep=''))
                par(mfrow=c(2,2))
                
                #- boxplot en fonction des categories
                categories <- unique(dat[,4])
                fc <- list() ; lfc <- list()
                for (category in categories) {{
                    idx <- which(dat[,4] == category)
                    fc[[category]] <- dat[idx,5]
                    lfc[[category]] <- log2(dat[idx,5])
                }}
                colors <- rep('lightgrey',length(categories))
                if ('Random' %in% categories) {{ colors[which(categories == 'Random')] <- 'darkblue' }}
                if ('PosEpromoter' %in% categories) {{ colors[which(categories == 'PosEpromoter')] <- 'darkgreen' }}
                if ('Positive' %in% categories) {{ colors[which(categories == 'Positive')] <- 'darkgreen' }}
                if ('Negative' %in% categories) {{ colors[which(categories == 'Negative')] <- 'darkred' }}
                boxplot(fc, pch=20, col=colors, main='Fold Change of categories', ylab='Fold Change',las=2)
                text(length(categories), 0.9*max(dat[,5]), labels=sprintf('Threshold :\n%3.2f',th_FoldChange), col='red')
                boxplot(lfc, pch=20, col=colors, main='Fold Change of categories', ylab='Fold Change [log2]', las=2)
                abline(h=log2(th_FoldChange), col='red', lty='dashed')

                
                #- ranked genomic regions based on their Fold Change
                plot(sort(dat[-which(dat[,4]=='Random' | dat[,4]=='Negative'),5]), pch=20, main='Activity of genomic regions\n(Random/Negative not included)', ylab='Fold Change', xlab='Ranked genomic regions')
                idx <- which(dat[,5] >= th_FoldChange & dat[,4] != 'Random' & dat[,4] != 'Negative')
                abline(v=nrow(dat)-length(idx), col='red', lty='dashed')
                text(nrow(dat)-length(idx), 0.9*max(dat[,5]), labels=length(idx), pos=2, offset=0, col='red')
                idx <- which(dat[,5] < th_FoldChange & dat[,4] != 'Random' & dat[,4] != 'Negative')
                text(length(idx)/2, 0.9*max(dat[,5]), labels=length(idx), col='black')
        
                #- info sur l'analyse
                plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
                text(x=0.5, y=1, 'INFORMATION', cex=1.5, pos=1, offset=0, col = 'black')
                text(x=0, y=0.70, '{wildcards.id_sample}', pos=4, offset=0, col = 'black')
                text(x=0, y=0.60, paste('Method: FDR=',fdr,sep=''), pos=4, offset=0, col = 'black')            
                dev.off()            
                
                dat_groups <- data.frame(dat[,1:4], groups)
                write.table(dat_groups, file=paste('{params.outdir}/Groups.FDR=',fdr,'.grp',sep=''), quote=F, row.names=F, col.names=F, sep='\t')
            }}
        }}
        
        
        ## CALCUL DU SEUIL EN FONCTION DU POINT D'INFLEXION
        
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
                abline(v= xPt,h= y_cutoff,lty=2,col=8)
                points(xPt,y_cutoff,pch=16,cex=0.9,col=2)
                abline(coef=c(b,slope),col=2)
                title(paste("x=",xPt,"\ny=",signif(y_cutoff,3),"\nFold over Median=",signif(y_cutoff/median(inputVector),3),"x\nFold over Mean=",signif(y_cutoff/mean(inputVector),3),"x",sep=""))
                axis(1,sum(inputVector==0),sum(inputVector==0),col.axis="pink",col="pink") #Number of regions with zero signal
            }}
            return(list(absolute=y_cutoff,overMedian=y_cutoff/median(inputVector),overMean=y_cutoff/mean(inputVector)))
        }}

        #this is an accessory function, that determines the number of points below a diagnoal passing through [x,yPt]
        numPts_below_line <- function(myVector,slope,x){{
            yPt <- myVector[x]
            b <- yPt-(slope*x)
            xPts <- 1:length(myVector)
            return(sum(myVector<=(xPts*slope+b)))
        }}
        #--
        #---  code from ROSE tool to dertermine super-enhancers  ---#

        # determination du seuil
        idx <- which(dat[,4] != 'Negative' & dat[,4] != 'Random')
        th_FoldChange <- calculate_cutoff(dat[idx,5], drawPlot=F)$absolute

        # identification des regions actives/inactives
        idx <- which(dat[,5] >= th_FoldChange)
        groups <- rep('Inactive',nrow(dat))
        groups[idx] <- 'Active'
        groups[which(is.na(dat[,5]))] <- 'NA'

        pdf('{output.pdf}')
        par(mfrow=c(2,2))

        #- boxplot en fonction des categories
        categories <- unique(dat[,4])
        fc <- list() ; lfc <- list()
        for (category in categories) {{
            idx <- which(dat[,4] == category)
            fc[[category]] <- dat[idx,5]
            lfc[[category]] <- log2(dat[idx,5])
        }}
        colors <- rep('lightgrey',length(categories))
        if ('Random' %in% categories) {{ colors[which(categories == 'Random')] <- 'darkblue' }}
        if ('PosEpromoter' %in% categories) {{ colors[which(categories == 'PosEpromoter')] <- 'darkgreen' }}
        if ('Positive' %in% categories) {{ colors[which(categories == 'Positive')] <- 'darkgreen' }}
        if ('Negative' %in% categories) {{ colors[which(categories == 'Negative')] <- 'darkred' }}
        boxplot(fc, pch=20, col=colors, main='Fold Change of categories', ylab='Fold Change', las=2)
        text(length(categories), 0.9*max(dat[,5]), labels=sprintf('Threshold :\n%3.2f',th_FoldChange), col='red')
        boxplot(lfc, pch=20, col=colors, main='Fold Change of categories', ylab='Fold Change [log2]', las=2)
        abline(h=log2(th_FoldChange), col='red', lty='dashed')
        
        #- ranked genomic regions based on their Fold Change
        plot(sort(dat[which(dat[,4]!='Random' & dat[,4]!='Negative'),5]), pch=20, main='Activity of genomic regions\n(Random/Negative not included)', ylab='Fold Change', xlab='Ranked genomic regions')

        idx <- which(dat[,5] >= th_FoldChange & dat[,4] != 'Random' & dat[,4] != 'Negative')
        abline(v=nrow(dat)-length(idx), col='red', lty='dashed')
        text(nrow(dat)-length(idx)/2, 0.9*max(dat[,5]), labels=length(idx), pos=2, offset=0, col='red')
        idx <- which(dat[,5] < th_FoldChange & dat[,4] != 'Random' & dat[,4] != 'Negative')
        text(length(idx)/2, 0.9*max(dat[,5]), labels=length(idx), col='black')
        
        print('debug_test4')
        
        #- info sur l'analyse
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
        text(x=0.5, y=1, 'INFORMATION', cex=1.5, pos=1, offset=0, col = 'black')
        text(x=0, y=0.70, '{wildcards.id_sample}', pos=4, offset=0, col = 'black')
        text(x=0, y=0.60, 'Method: Inflexion point', pos=4, offset=0, col = 'black')
        dev.off()            

        dat <- data.frame(dat[,1:4], groups)
        write.table(dat, file='{output.groups}', quote=F, row.names=F, col.names=F, sep='\t')
        """)

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
    params:
    run:
        shell("cp {input.pdf} {output.pdf}")
        R("""
        bed_crm <- read.table('{input.bed_crm}', stringsAsFactors=F, sep="\t")
        fpkm <- read.table('{input.fpkm}', stringsAsFactors=F)
        fpkm_input <- read.table('{input.fpkm_input}', stringsAsFactors=F)
        fc <- read.table('{input.fc}', stringsAsFactors=F)
        group <- read.table('{input.group}', stringsAsFactors=F)
     
        # Because some of my input bed are bed4 and some are bed6.
        if (dim(bed_crm)[2] == 4){{
            bed_crm[,5] <- '.'
            bed_crm[,6] <- 'NA'
            }}

        dat <- data.frame(fpkm, fpkm_input[,5], fc[,5], group[,5], bed_crm[,6])
        #group_labels <- sapply(strsplit(group_files,'\\\.'), function(x){{ nb=length(x) ; return(paste('group_',x[nb-2],'.',x[nb-1],sep='')) }})
        colnames(dat) <- c('chr','start','end','category','fpkm','fpkm_input','fold_change','group_inflexionPoint','gene_name')
        
        #if (file.exists('data/CRMs/CRMs_geneAssociation.bed_notOrdered') == TRUE) {{
        #    genes <- read.table('data/CRMs/CRMs_geneAssociation.bed_notOrdered', stringsAsFactors=F)
        #    str_fpkm <- apply(fpkm[,1:3],1,paste,collapse='_')
        #    str_genes <- apply(genes[,1:3],1,paste,collapse='_')
        #    idx <- match(str_fpkm, str_genes)
        #    genes <- genes[idx,]
        #    dat <- data.frame(dat, genes=genes[,4], location=genes[,5])
        #}}

        write.table(dat, file='{output.tsv}', quote=F, row.names=F, col.names=T, sep='\t')
        """)

#localrules: targets, ln_seq_id_to_sample_id_bam

#rule target:
#    threads: 1
#    message: "-- Rule target completed. --"
#    input:
#        # bamCoverage ext314 (bigiwig)
#        ## with removed duplicates:
#        expand("out/10_314/7/5_mm9/3/2/1/{id}.bw",
#            id=[
#                "run155_CapSTAR-seq_mSilencers_mT_DHS_PGK_P5424_rep1",
#                "run155_CapSTAR-seq_mSilencers_mT_DHS_PGK_input",
#                "run170_CapSTAR-seq_mSilencers_mT_DHS_PGK_P5424_rep2",
#                "run170_CapSTAR-seq_mT_DHS_EFIA_P5424_rep1",
#                "run170_CapSTAR-seq_mT_DHS_EFIA_P5424_rep2",
#                "run170_CapSTAR-seq_mT_DHS_EFIA_P5424_rep3",
#                "run170_CapSTAR-seq_mT_DHS_EFIA_input"]),
#        # FoldChange FPKM on CRMs
#        ## with removed duplicates:
#        "out/19_mTDHS/14_314/13/7/5_mm9/3/2/1/run155_CapSTAR-seq_mSilencers_mT_DHS_PGK_P5424_rep1_over_run155_CapSTAR-seq_mSilencers_mT_DHS_PGK_input.allData.tsv",
#        "out/19_mTDHS/14_314/13/7/5_mm9/3/2/1/run170_CapSTAR-seq_mSilencers_mT_DHS_PGK_P5424_rep2_over_run155_CapSTAR-seq_mSilencers_mT_DHS_PGK_input.allData.tsv",
#        "out/19_mTDHS/14_314/13/7/5_mm9/3/2/1/run170_CapSTAR-seq_mT_DHS_EFIA_P5424_rep1_over_run170_CapSTAR-seq_mT_DHS_EFIA_input.allData.tsv",
#        "out/19_mTDHS/14_314/13/7/5_mm9/3/2/1/run170_CapSTAR-seq_mT_DHS_EFIA_P5424_rep2_over_run170_CapSTAR-seq_mT_DHS_EFIA_input.allData.tsv",
#        "out/19_mTDHS/14_314/13/7/5_mm9/3/2/1/run170_CapSTAR-seq_mT_DHS_EFIA_P5424_rep3_over_run170_CapSTAR-seq_mT_DHS_EFIA_input.allData.tsv"
#    shell:"""
#    WDIR=`pwd`
#    mkdir -p result/capstarseq_2017_01_05
#    cd result/capstarseq_2017_01_05
#    ln -s ../../out/10_314/7/5_mm9/3/2/1/ bw
#    ln -s ../../out/19_mTDHS/14_314/13/7/5_mm9/3/2/1/ tab
#    #echo "workflow completed" | mail -s "workflow completed" guillaume.charbonnier@outlook.com
#    #chown -R Charbonnier:sacapus .
#    #chmod -R g+rwX .
#    """


##def seq_id_to_sample_id(wildcards):
##	"""
##	Created: 2016-05-27
##	"""
##	sample = wildcards['sample']
##	exp = wildcards['exp']
##	d={}
##	for row in csv.DictReader(open('code/snakemake/samples.tsv'),delimiter='\t'):
##		d[row['Sample_name']] = row['File_name']
##	return "inp/fastq/" + exp + "/" + d[sample]
#
#def seq_id_to_sample_id_for_ln_bam(wildcards):
#    """
#    Created: 2016-05-27
#    """
#    sample = wildcards['sample']
#    #print(sample)
#    exp = wildcards['exp']
#    index = wildcards['index']
#    ext = wildcards['ext']
#    #print(ext)
#    d={}
#    for row in csv.DictReader(open('code/snakemake/samples.tsv'),delimiter='\t'):
#        d[row['Sample_name']] = row['File_name']
#    return "out/samtools/sam_to_bam/" + index + "/" + exp + "/" + d[sample] + "." + ext
#
#

# This r14 is translated into awk_extend_reads in awk.rules.
#rule r14_awk_extReads:
#    """
#    Taken from Aurelien's workflow.
#    Added double escape for tabulation in awk.
#    Added check condition if chromosomes are like "chr3" or only "3".
#    """
#    input:  bed="out/{id}.bed"
#    output: bed="out/14_{extReads}/{id}.bed"
#    shell:"""
#    size_fragment={wildcards.extReads}
#    awk -v FRAG_SIZE=$size_fragment 'BEGIN{{ OFS="\\t" }}{{
#        if ($1 ~ "^chr*"){{
#            if ($6 == "+"){{ $3 = $2 + FRAG_SIZE ; print $0 }}
#            else if ($6 == "-"){{ $2 = $3 - FRAG_SIZE ; if ($2 > 0) {{ print $0 }} }}
#            }}
#        else {{
#            if ($6 == "+"){{ $3 = $2 + FRAG_SIZE ; print "chr"$0 }}
#            else if ($6 == "-"){{ $2 = $3 - FRAG_SIZE ; if ($2 > 0) {{ print "chr"$0 }} }}
#            }}
#        }}' {input.bed} > {output.bed}
#    """
#

#
#rule capstarrseq_bedtools_coverage_fpkm_input:
#    """
#    Test:
#        out/15_mTDHS/awk/extend_reads_314/bedtools/bamtobed/samtools/sam_to_bam/bowtie2/pe_mm9/sickle/pe_-t_sanger_-q_20/gunzip/ln/rename_run178_tgml/BAC_mT_DHS_PGK_p5424_rep1.filtered_CRMs.bed
#    """
#    input:
#        bedtools="opt/miniconda/envs/bedtools/bin/bedtools",
#        bed_reads="out/{id_bam_to_bed}/{id}.bed",
#        bed_crms="inp/bed/crms/{crm_type}.bed",
#        flagstat="out/{id}.flagstat.txt"
#    output:
#        tsv_cov_unfiltered="out/15_{crm_type}/{id_bam_to_bed}/{id}.coverage_unfiltered.tsv",
#        tsv_fpkm="out/15_{crm_type}/{id_bam_to_bed}/{id}.FPKM.tsv",
#        tsv_fpkm_unfiltered="out/15_{crm_type}/{id_bam_to_bed}/{id}.FPKM_unfiltered.tsv",
#        bed_crms="out/15_{crm_type}/{id_bam_to_bed}/{id}.filtered_CRMs.bed"
#    params:
#        fpkm_threshold='0'
#    wildcard_constraints:
#        id_bam_to_bed="bedtools/bamtobed|awk/extend_reads_[0-9]+/bedtools/bamtobed",
#        crm_type="mTDHS|hProm|hProm_posEprom"
#    run:
#        shell("""# 1. Coverage
#        {input.bedtools} coverage -a {input.bed_crms} -b {input.bed_reads} > {output.tsv_cov_unfiltered}
#        #cp {output.tsv_cov_unfiltered} tmp/debug_tsv_cov_unfiltered.tsv
#        """)
#        shell("""# 2. Coverage to FPKM
#        nb_mReads=`grep "mapped (" {input.flagstat} | awk '{{print $1}}'`
#        awk -v NB_READS=$nb_mReads '{{ fpkm = ($(NF-3) / (($(NF-1))/1000)) / (NB_READS/1000000) ; print $1"\\t"$2"\\t"$3"\\t"$4"\\t"fpkm }}' {output.tsv_cov_unfiltered} > {output.tsv_fpkm_unfiltered}
#        #cp {output.tsv_fpkm_unfiltered} tmp/debug_capstarrseq_bedtools_coverage_fpkm_input.tsv
#        """)
#        R("""
#        # 3. CRM filtering based on FPKM threshold (input FPKM < 1)
#        fpkm_all <- read.table('{output.tsv_fpkm_unfiltered}', stringsAsFactors=F)
#        idx <- which(fpkm_all[,5] >= {params.fpkm_threshold})
#        fpkm_all <- fpkm_all[idx,]
#        write.table(fpkm_all, file='{output.tsv_fpkm}', quote=F, row.names=F, col.names=F, sep='\t')
#
#        # 4. Producing filtered CRM file
#        crms_all <- read.table('{input.bed_crms}', stringsAsFactors=F)
#        #fpkm_all <- read.table('{output.tsv_fpkm}', stringsAsFactors=F)
#        str_crms_all <- apply(crms_all[,1:3], 1, paste, collapse='_')
#        str_fpkm_all <- apply(fpkm_all[,1:3], 1, paste, collapse='_')
#        idx <- match(str_fpkm_all, str_crms_all)
#        write.table(crms_all[idx,], file='{output.bed_crms}', quote=F, row.names=F, col.names=F, sep='\t')
#        """)
#
#
##def input_crms_bedtools_coverage_fpkm_sample(wildcards):
##    crm_type = wildcards['crm_type']
##    index = wildcards['index']
##    extReads = wildcards['extReads']
##    exp = wildcards['exp']
##
##    if exp == 'run149_150_CapSTAR-seq_hPromoters':
##        exp_sample_input = "run107_CapSTAR-seq_hpromoter/Input"
##        
##    elif exp == 'run107_CapSTAR-seq_hpromoter':
##        exp_sample_input = "run107_CapSTAR-seq_hpromoter/Input"
##        
##    elif exp == 'run155_CapSTAR-seq_mSilencers':
##        exp_sample_input = "run155_CapSTAR-seq_mSilencers/mT_DHS_PGK_input"
##
##    bed_crms="out/bedtools/coverage/fpkm_input/CRMs_"+crm_type+"/extReads"+extReads+"/"+index+"/"+exp_sample_input+"/filtered_CRMs.bed"
##    return bed_crms
#
#rule r16_bedtools_coverage_fpkm_sample:
#    input:
#        bedtools="opt/miniconda/envs/bedtools/bin/bedtools",
#        bed_reads="out/{id_bam_to_bed}/{id}/{id_sample}.bed",
#        bed_crms="out/15_{crm_type}/{id_bam_to_bed}/{id}/{id_input}.filtered_CRMs.bed",
#        flagstat="out/{id}/{id_sample}.flagstat.txt"
#    output:
#        tsv_cov="out/16_{crm_type}/{id_bam_to_bed}/{id}/{id_sample}_on_{id_input}_filtered_CRMs.coverage.tsv",
#        tsv_fpkm="out/16_{crm_type}/{id_bam_to_bed}/{id}/{id_sample}_on_{id_input}_filtered_CRMs.FPKM.tsv"
#    wildcard_constraints:
#        id_bam_to_bed="bedtools/bamtobed|awk/extend_reads_[0-9]+/bedtools/bamtobed",
#        crm_type="mTDHS|hProm|hProm_posEprom"
#    shell:"""
#    # 1. Coverage
#    {input.bedtools} coverage -a {input.bed_crms} -b {input.bed_reads} > {output.tsv_cov}
#
#    # 2. Coverage to FPKM
#    nb_mReads=`grep "mapped (" {input.flagstat} | awk '{{print $1}}'`
#    awk -v NB_READS=$nb_mReads '{{ fpkm = ($(NF-3) / (($(NF-1))/1000)) / (NB_READS/1000000) ; print $1"\\t"$2"\\t"$3"\\t"$4"\\t"fpkm }}' {output.tsv_cov} > {output.tsv_fpkm}
#    """
#
###def input_bed_crms_bedtools_coverage(wildcards):
###    ori_filt = wildcards['ori_filt']
###    crm_type = wildcards['crm_type']
###    index = wildcards['index']
###    extReads = wildcards['extReads']
###    exp = wildcards['exp']
###
###    if ori_filt == "ori":
###        bed_crms="inp/bed/crms/"+crm_type+".bed"
###
###    elif ori_filt == "filt":
###        if exp == 'run149_150_CapSTAR-seq_hPromoters':
###            exp_sample_input = "run107_CapSTAR-seq_hpromoter/Input"
###        
###        elif exp == 'run107_CapSTAR-seq_hpromoter':
###            exp_sample_input = "run107_CapSTAR-seq_hpromoter/Input"
###        
###        elif exp == 'run155_CapSTAR-seq_mSilencers':
###            exp_sample_input = "run155_CapSTAR-seq_mSilencers/mT_DHS_PGK_input"
###
###        bed_crms="out/r/crm_filtering/"+crm_type+"/extReads"+extReads+"/"+index+"/"+exp_sample_input+".bed"
###    print('print bed_crms:')
###    print(bed_crms)
###    return bed_crms
###
###
###rule bedtools_coverage:
###    """
###    Note: expand CRMs_all later.
###    "out/bedtools/coverage/CRMs_hProm_filt/extReads314/hg19/run149_150_CapSTAR-seq_hPromoters/RPMI_starr_hP_rep1.tsv",$
###    """
###    input:
###        bedtools="opt/miniconda/envs/py35/bin/bedtools",
###        bed_reads="out/awk/extReads{extReads}/{index}/{exp}/{id}.bed",
###        bed_crms=input_bed_crms_bedtools_coverage
###    output: tsv="out/bedtools/coverage/CRMs_{crm_type}_{ori_filt}/extReads{extReads}/{index}/{exp}/{id}.tsv"
###    wildcard_constraints:
###        crm_type="hProm|IGMM|mTDHS",
###        ori_filt="ori|filt",
###        extReads="[0-9]+",
###        index="hg19|mm9"
###    shell:"""
###    {input.bedtools} coverage -a {input.bed_crms} -b {input.bed_reads} > {output.tsv}
###    """
###
###rule awk_compute_fpkm:
###    """
###    Taken from Aurelien's workflow
###    """
###    input:
###        txt="out/samtools/flagstat/{index}/{exp}/{id}.txt",
###        tsv="out/bedtools/coverage/CRMs_{crm_type}_{ori_filt}/extReads{extReads}/{index}/{exp}/{id}.tsv"
###    output: tsv="out/awk/compute_fpkm/CRMs_{crm_type}_{ori_filt}/extReads{extReads}/{index}/{exp}/{id}.tsv"
###    wildcard_constraints: ori_filt = "ori|filt"
###    shell:"""
###    nb_mReads=`grep "mapped (" {input.txt} | awk '{{print $1}}'`
###    #nb_mReads=`awk 'BEGIN{{nb_mReads=0}} {{nb_mReads=nb_mReads+$5}} END{{print nb_mReads}}' #params.outdir#/analysis/coverage`
###    
###    awk -v NB_READS=$nb_mReads '{{ fpkm = ($(NF-3) / (($(NF-1))/1000)) / (NB_READS/1000000) ; print $1"\\t"$2"\\t"$3"\\t"$4"\\t"fpkm }}' {input.tsv} > {output.tsv}
###    """
###
###
###def input_r_crm_filtering(wildcards):
###    crm_type = wildcards['crm_type']
###    extReads = wildcards['extReads']
###    index = wildcards['index']
###    exp_sample_input = wildcards["exp_sample_input"]
###
###    if crm_type == 'hProm':
###        exp_sample_input = "run107_CapSTAR-seq_hpromoter/Input"
###    
###    elif crm_type == 'mTDHS':
###        exp_sample_input = "run155_CapSTAR-seq_mSilencers/mT_DHS_PGK_input"
###    
###    #elif crm_type == 'IGMM':
###    #    exp_sample_input = "run155_CapSTAR-seq_mSilencers/mT_DHS_PGK_in
###    
###    tsv_path = ''.join(["out/awk/compute_fpkm/CRMs_", crm_type, "_ori/extReads", extReads, "/", index, "/", exp_sample_input, ".tsv"])
###    return tsv_path
###
###
###rule r_crm_filtering:
###    """Filter input CRMs based on FPKM threshold observed in input sample"""
###    input:
###        tsv=input_r_crm_filtering,
###        bed_crm="inp/bed/crms/{crm_type}.bed"
###    output:
###        bed_crm="out/r/crm_filtering/{crm_type}/{extReads}/{index}/{exp_sample_input}.bed"
###    params: input_fpkm_threshold='1'
###    run:
###        R("""
###        #- filtrage des CRMs selon le FPKM (input FPKM < 1)
###        fpkm_all <- read.table('{input.tsv}', stringsAsFactors=F)
###        idx <- which(fpkm_all[,5] >= {params.input_fpkm_threshold})
###        fpkm_all <- fpkm_all[idx,]
###        #write.table(fpkm_all, file='{params.outdir}/analysis/FPKM_all.fpkm', quote=F, row.names=F, col.names=F, sep='\t')
###
###        #- elimination des CRMs filtres dans les listes des CRMs etudies
###        crms_all <- read.table('{input.bed_crm}', stringsAsFactors=F)
###        #fpkm_all <- read.table('{params.outdir}/analysis/FPKM_all.fpkm', stringsAsFactors=F)
###        str_crms_all <- apply(crms_all[,1:3], 1, paste, collapse='_')
###        str_fpkm_all <- apply(fpkm_all[,1:3], 1, paste, collapse='_')
###        idx <- match(str_fpkm_all, str_crms_all)
###        write.table(crms_all[idx,], file='{output.bed_crm}', quote=F, row.names=F, col.names=F, sep='\t')
###        """)
###
###
##def input_r_fold_change(wildcards):
##    exp = wildcards['exp']
##    crm_type = wildcards['crm_type']
##    extReads = wildcards['extReads']
##    index = wildcards['index']
##
##    if exp == 'run149_150_CapSTAR-seq_hPromoters':
##        exp_sample_input = "run107_CapSTAR-seq_hpromoter/Input"
##    
##    elif exp == 'run107_CapSTAR-seq_hpromoter':
##        exp_sample_input = "run107_CapSTAR-seq_hpromoter/Input"
##
##    elif exp == 'run155_CapSTAR-seq_mSilencers':
##        exp_sample_input = "run155_CapSTAR-seq_mSilencers/mT_DHS_PGK_input"
##
##    else:
##        print('add condition to input_r_fold_change function')
##
##    path_input = ''.join(["out/bedtools/coverage/fpkm_input/CRMs_", crm_type, "/extReads", extReads, "/", index, "/", exp_sample_input, "/FPKM.tsv"])
##    return path_input
#
#
## "out/17/14_314/.../rep1_over_rep2.foldChange.tsv"
#rule r17_fold_change:
#    """
#    Calcul des Fold Change dans les echantillons
#    """
#    input:
#        tsv_fpkm_sample = "out/16_{crm_type}/{id}/{id_sample}_on_{id_input}_filtered_CRMs.FPKM.tsv",
#        tsv_fpkm_input = "out/15_{crm_type}/{id}/{id_input}.FPKM.tsv"
#    output:
#        fc = "out/17_{crm_type}/{id}/{id_sample}_over_{id_input}.foldChange.tsv"
#    wildcard_constraints:
#        #id_bam_to_bed="bedtools/bamtobed|awk/extend_reads_[0-9]+/bedtools/bamtobed",
#        crm_type="mTDHS|hProm|hProm_posEprom"
#    run:
#        R("""
#        fpkm_input <- read.table('{input.tsv_fpkm_input}', stringsAsFactors=F)
#        fpkm_sample <- read.table('{input.tsv_fpkm_sample}', stringsAsFactors=F)
#        fc <- fpkm_sample[,5] / fpkm_input[,5]
#        dat <- data.frame(fpkm_sample[,1:4], fc)
#        write.table(dat, file='{output.fc}', quote=F, row.names=F, col.names=F, sep='\t')
#        """)
#
#rule r18_grouping_crms:
#    """
#    Creation des groupes de CRMs (inactifs, actifs) sur la base des Fold Change et des categories
#    """
#    input:
#        fc = "out/{id}/{id_sample}_over_{id_input}.foldChange.tsv"
#    output:
#        groups = "out/18/{id}/{id_sample}_over_{id_input}.inflexionPointGroups.tsv",
#        pdf = "out/18/{id}/{id_sample}_over_{id_input}.inflexionPointGroups.pdf"
#    params:
#        outdir = "out/18/{id}/",
#        th_negative_1 = "0.05",
#        th_negative_2 = "0.01"
#    run:
#        R("""
#        dat <- read.table('{input.fc}', stringsAsFactors=F)
#    
#        ## CALCUL DU SEUIL EN FONCTION D'UN FDR
#        for (fdr in c({params.th_negative_1},{params.th_negative_2})) {{
#            # Determination du seuil
#            idx <- which(dat[,4] == 'Negative')
#            # Skipping FDR computation if no Negative set is provided.
#            if (length(idx)!=0) {{
#                P <- ecdf(dat[idx,5])
#                th_FoldChange <- quantile(P,probs=1-fdr)
#
#                # identification des regions actives/inactives
#                idx <- which(dat[,5] >= th_FoldChange)
#                groups <- rep('Inactive',nrow(dat))
#                groups[idx] <- 'Active'
#                groups[which(is.na(dat[,5]))] <- 'NA'
#                
#                pdf(paste('{params.outdir}/Groups.FDR=',fdr,'.pdf',sep=''))
#                par(mfrow=c(2,2))
#                
#                #- boxplot en fonction des categories
#                categories <- unique(dat[,4])
#                fc <- list() ; lfc <- list()
#                for (category in categories) {{
#                    idx <- which(dat[,4] == category)
#                    fc[[category]] <- dat[idx,5]
#                    lfc[[category]] <- log2(dat[idx,5])
#                }}
#                colors <- rep('lightgrey',length(categories))
#                if ('Random' %in% categories) {{ colors[which(categories == 'Random')] <- 'darkblue' }}
#                if ('Positive' %in% categories) {{ colors[which(categories == 'Positive')] <- 'darkgreen' }}
#                if ('Negative' %in% categories) {{ colors[which(categories == 'Negative')] <- 'darkred' }}
#                boxplot(fc, pch=20, col=colors, main='Fold Change of categories', ylab='Fold Change')
#                text(length(categories), 0.9*max(dat[,5]), labels=sprintf('Threshold :\n%3.2f',th_FoldChange), col='red')
#                boxplot(lfc, pch=20, col=colors, main='Fold Change of categories', ylab='Fold Change [log2]')
#                abline(h=log2(th_FoldChange), col='red', lty='dashed')
#
#                
#                #- ranked genomic regions based on their Fold Change
#                plot(sort(dat[-which(dat[,4]=='Random' | dat[,4]=='Negative'),5]), pch=20, main='Activity of genomic regions\n(Random/Negative not included)', ylab='Fold Change', xlab='Ranked genomic regions')
#                idx <- which(dat[,5] >= th_FoldChange & dat[,4] != 'Random' & dat[,4] != 'Negative')
#                abline(v=nrow(dat)-length(idx), col='red', lty='dashed')
#                text(nrow(dat)-length(idx), 0.9*max(dat[,5]), labels=length(idx), pos=2, offset=0, col='red')
#                idx <- which(dat[,5] < th_FoldChange & dat[,4] != 'Random' & dat[,4] != 'Negative')
#                text(length(idx)/2, 0.9*max(dat[,5]), labels=length(idx), col='black')
#        
#                #- info sur l'analyse
#                plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
#                text(x=0.5, y=1, 'INFORMATION', cex=1.5, pos=1, offset=0, col = 'black')
#                text(x=0, y=0.70, '{wildcards.id_sample}', pos=4, offset=0, col = 'black')
#                text(x=0, y=0.60, paste('Method: FDR=',fdr,sep=''), pos=4, offset=0, col = 'black')            
#                dev.off()            
#                
#                dat_groups <- data.frame(dat[,1:4], groups)
#                write.table(dat_groups, file=paste('{params.outdir}/Groups.FDR=',fdr,'.grp',sep=''), quote=F, row.names=F, col.names=F, sep='\t')
#            }}
#        }}
#        
#        
#        ## CALCUL DU SEUIL EN FONCTION DU POINT D'INFLEXION
#        
#        #---  code from ROSE tool to dertermine super-enhancers  ---#
#        #--
#        #This function calculates the cutoff by sliding a diagonal line and finding where it is tangential (or as close as possible)
#        calculate_cutoff <- function(inputVector, drawPlot=TRUE,...){{
#            inputVector <- sort(inputVector)
#            inputVector[inputVector<0]<-0 #set those regions with more control than ranking equal to zero
#            slope <- (max(inputVector)-min(inputVector))/length(inputVector) #This is the slope of the line we want to slide. This is the diagonal.
#            xPt <- floor(optimize(numPts_below_line,lower=1,upper=length(inputVector),myVector= inputVector,slope=slope)$minimum) #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
#            y_cutoff <- inputVector[xPt] #The y-value at this x point. This is our cutoff.
#    
#            if(drawPlot){{  #if TRUE, draw the plot
#                plot(1:length(inputVector), inputVector,type="l",...)
#                b <- y_cutoff-(slope* xPt)
#                abline(v= xPt,h= y_cutoff,lty=2,col=8)
#                points(xPt,y_cutoff,pch=16,cex=0.9,col=2)
#                abline(coef=c(b,slope),col=2)
#                title(paste("x=",xPt,"\ny=",signif(y_cutoff,3),"\nFold over Median=",signif(y_cutoff/median(inputVector),3),"x\nFold over Mean=",signif(y_cutoff/mean(inputVector),3),"x",sep=""))
#                axis(1,sum(inputVector==0),sum(inputVector==0),col.axis="pink",col="pink") #Number of regions with zero signal
#            }}
#            return(list(absolute=y_cutoff,overMedian=y_cutoff/median(inputVector),overMean=y_cutoff/mean(inputVector)))
#        }}
#
#        #this is an accessory function, that determines the number of points below a diagnoal passing through [x,yPt]
#        numPts_below_line <- function(myVector,slope,x){{
#            yPt <- myVector[x]
#            b <- yPt-(slope*x)
#            xPts <- 1:length(myVector)
#            return(sum(myVector<=(xPts*slope+b)))
#        }}
#        #--
#        #---  code from ROSE tool to dertermine super-enhancers  ---#
#
#        # determination du seuil
#        idx <- which(dat[,4] != 'Negative' & dat[,4] != 'Random')
#        th_FoldChange <- calculate_cutoff(dat[idx,5], drawPlot=F)$absolute
#
#        # identification des regions actives/inactives
#        idx <- which(dat[,5] >= th_FoldChange)
#        groups <- rep('Inactive',nrow(dat))
#        groups[idx] <- 'Active'
#        groups[which(is.na(dat[,5]))] <- 'NA'
#
#        pdf('{output.pdf}')
#        par(mfrow=c(2,2))
#
#        #- boxplot en fonction des categories
#        categories <- unique(dat[,4])
#        fc <- list() ; lfc <- list()
#        for (category in categories) {{
#            idx <- which(dat[,4] == category)
#            fc[[category]] <- dat[idx,5]
#            lfc[[category]] <- log2(dat[idx,5])
#        }}
#        colors <- rep('lightgrey',length(categories))
#        if ('Random' %in% categories) {{ colors[which(categories == 'Random')] <- 'darkblue' }}
#        if ('Positive' %in% categories) {{ colors[which(categories == 'Positive')] <- 'darkgreen' }}
#        if ('Negative' %in% categories) {{ colors[which(categories == 'Negative')] <- 'darkred' }}
#        boxplot(fc, pch=20, col=colors, main='Fold Change of categories', ylab='Fold Change')
#        text(length(categories), 0.9*max(dat[,5]), labels=sprintf('Threshold :\n%3.2f',th_FoldChange), col='red')
#        boxplot(lfc, pch=20, col=colors, main='Fold Change of categories', ylab='Fold Change [log2]')
#        abline(h=log2(th_FoldChange), col='red', lty='dashed')
#
#        #- ranked genomic regions based on their Fold Change
#        plot(sort(dat[-which(dat[,4]=='Random' | dat[,4]=='Negative'),5]), pch=20, main='Activity of genomic regions\n(Random/Negative not included)', ylab='Fold Change', xlab='Ranked genomic regions')
#        idx <- which(dat[,5] >= th_FoldChange & dat[,4] != 'Random' & dat[,4] != 'Negative')
#        abline(v=nrow(dat)-length(idx), col='red', lty='dashed')
#        text(nrow(dat)-length(idx)/2, 0.9*max(dat[,5]), labels=length(idx), pos=2, offset=0, col='red')
#        idx <- which(dat[,5] < th_FoldChange & dat[,4] != 'Random' & dat[,4] != 'Negative')
#        text(length(idx)/2, 0.9*max(dat[,5]), labels=length(idx), col='black')
#        
#        #- info sur l'analyse
#        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
#        text(x=0.5, y=1, 'INFORMATION', cex=1.5, pos=1, offset=0, col = 'black')
#        text(x=0, y=0.70, '{wildcards.id_sample}', pos=4, offset=0, col = 'black')
#        text(x=0, y=0.60, 'Method: Inflexion point', pos=4, offset=0, col = 'black')
#        dev.off()            
#
#        dat <- data.frame(dat[,1:4], groups)
#        write.table(dat, file='{output.groups}', quote=F, row.names=F, col.names=F, sep='\t')
#        """)
#
#rule r19_merge_all_data:
#    input:
#        fpkm = "out/16_{crm_type}/{id}/{id_sample}_on_{id_input}_filtered_CRMs.FPKM.tsv",
#        fpkm_input = "out/15_{crm_type}/{id}/{id_input}.FPKM.tsv",
#        fc = "out/17_{crm_type}/{id}/{id_sample}_over_{id_input}.foldChange.tsv",
#        group = "out/18/17_{crm_type}/{id}/{id_sample}_over_{id_input}.inflexionPointGroups.tsv"
#    output:
#        "out/19_{crm_type}/{id}/{id_sample}_over_{id_input}.allData.tsv"
#    params:
#    run:
#        R("""
#        fpkm <- read.table('{input.fpkm}', stringsAsFactors=F)
#        fpkm_input <- read.table('{input.fpkm_input}', stringsAsFactors=F)
#        fc <- read.table('{input.fc}', stringsAsFactors=F)
#        group <- read.table('{input.group}', stringsAsFactors=F)
#       
#        dat <- data.frame(fpkm, fpkm_input[,5], fc[,5], group[,5])
#        #group_labels <- sapply(strsplit(group_files,'\\\.'), function(x){{ nb=length(x) ; return(paste('group_',x[nb-2],'.',x[nb-1],sep='')) }})
#        colnames(dat)[1:8] <- c('chr','start','end','category','fpkm','fpkm_input','fold_change','group_inflexionPoint')
#        
#        #if (file.exists('data/CRMs/CRMs_geneAssociation.bed_notOrdered') == TRUE) {{
#        #    genes <- read.table('data/CRMs/CRMs_geneAssociation.bed_notOrdered', stringsAsFactors=F)
#        #    str_fpkm <- apply(fpkm[,1:3],1,paste,collapse='_')
#        #    str_genes <- apply(genes[,1:3],1,paste,collapse='_')
#        #    idx <- match(str_fpkm, str_genes)
#        #    genes <- genes[idx,]
#        #    dat <- data.frame(dat, genes=genes[,4], location=genes[,5])
#        #}}
#
#        write.table(dat, file='{output}', quote=F, row.names=F, col.names=T, sep='\t')
#        """)
#
#
#rule fastQC:
#    """
#    Run fastQC control analysis on trimmed data.
#    """
#    input:
#        fastq="out/{source}/{exp_sample}.fastq",
#        fastqc="opt/FastQC/fastqc"
#    output: report="result/fastqc/{source, sickle|fastq}/{exp_sample}/"
#    threads: 1
#    shell:"""
#    {input.fastqc} --quiet --outdir {output.report} -f fastq {input.fastq}
#    """
#
#rule fastQC_bam:
#    """
#    Run fastQC control analysis on aligned data.
#    source="samtools/sam_to_bam/{index}/"
#    """
#    input:
#        bam="out/{source}/{index}/{exp_sample}.bam",
#        fastqc="opt/FastQC/fastqc"
#    output: report="result/fastqc/{source, samtools/sam_to_bam|picard/MarkDuplicates}/{index}/{exp_sample}/"
#    threads: 1
#    shell:"""
#    mkdir --parents {output.report}
#    {input.fastqc} --quiet --outdir {output.report} -f bam {input.bam}
#    """
#
#
#rule paired_end_to_single_end:
#        """
#        Created: 2016-07-18 18h13
#        Select only the reads that are first in a pair.
#        Bam produced contain reads and consistent size but nothing is displayed in IGV. I am not sure macs1.4 will take them as single-end for peak-calling so I will compare with macs2 paired-end peak-calling.
#        """
#        input:  samtools="opt/samtools-1.3.1/samtools",\
#                bam="out/ln/bam/{index}/{run}/{sample}.bam"
#        output: bam="out/samtools/pe_to_se/{index}/{run}/{sample}.bam"
#        threads: 2
#        shell:"""
#        {input.samtools} view -bSh -@ {threads} -f 0x0040 {input.bam} | \
#        {input.samtools} sort -@ {threads} -m 10G - > {output.bam}
#        {input.samtools} index {output.bam}
#        """
#
#rule paired_end_to_single_end_alter_sam_flag:
#        """
#        Created: 2016-07-29 16h30
#        Select only the reads that are first in a pair.
#        Bam produced contain reads and consistent size but nothing is displayed in IGV. I am not sure macs1.4 will take them as single-end for peak-calling so I will compare with macs2 paired-end peak-calling.
#        https://broadinstitute.github.io/picard/explain-flags.html
#        Try to substitute flags:
#69 69
#90 90
#113 16
#65 0
#73 0
#77 77
#81 16
#83 16
#89 16
#97 0
#99 0
#
#samtools view -h S001352_Input_Nut_wt-15015.bam | gawk 'BEGIN{FS=OFS="\t"}{if ($0 ~ /^@/) print; else {f2=and($2,16);l=$1"\t"f2;for (i=3;i<=NF;i++)l=l"\t"$i;print l}}' 
#        """
#        input:  samtools="opt/samtools-1.3.1/samtools",\
#                bam="out/ln/bam/{index}/{run}/{sample}.bam"
#        output: bam="out/samtools/pe_to_se_alter_sam_flag/{index}/{run}/{sample}.bam"
#        threads: 1
#        shell:"""
#        {input.samtools} view -h -@ {threads} {input.bam} |\
#        awk 'BEGIN{{FS=OFS="\\t"}}{{
#             # Printing the header
#             if ($0 ~ /^@/) print;
#             # Selecting lines that are first in pair
#             else if(and($2,64)==64){{
#                  # if read is reverse strand, then its bitwise flag is changed to reverse strand only (16), else forward strand only (0).
#                  f2=and($2,16)
#
#                  l=$1"\\t"f2;for (i=3;i<=NF;i++)l=l"\\t"$i;
#                  print l
#                  }}
#             }}' >tmp_awk_edit.sam #|\
#        #{input.samtools} view -bh - |\
#        #{input.samtools} sort -@ {threads} -m 10G - > {output.bam}
#        #{input.samtools} index {output.bam}
#        """
#
#
##def seq_id_to_sample_id_for_ln_bam_pe_to_se(wildcards):
##	"""
##	Created: 2016-05-27
##	"""
##	sample = wildcards['sample']
##	print(sample)
##	exp = wildcards['exp']
##	index = wildcards['index']
##	ext = wildcards['ext']
##	print(ext)
##	d={}
##	for row in csv.DictReader(open('code/snakemake/samples.tsv'),delimiter='\t'):
##		d[row['Sample_name']] = row['File_name']
##	return "out/samtools/pe_to_se/" + index + "/" + exp + "/" + d[sample] + "." + ext
##
##
##rule ln_seq_id_to_sample_id_bam:
##	input:	 seq_id_to_sample_id_for_ln_bam_pe_to_se
##	output: "out/ln/bam/{index}/{exp}/{sample, [a-zA-Z0-9_-]+}_pe_to_se.{ext, bam|bam.bai}"
##	shell:"""
##	ln {input} {output}
##	"""
##
##rule extract_and_convert_input_seq_id_to_explicit_sample_id:
##	input: seq_id_to_sample_id
##	output: "out/fastq/{exp}/{sample}.fastq"
##	shell:"""
##	gunzip --stdout {input} > {output}
##	"""
#
