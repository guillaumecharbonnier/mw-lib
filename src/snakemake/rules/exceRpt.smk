#!/bin/bash
#source activate /gpfs/tgml/nin/anaconda3/envs/exceRpt/
#conda activate exceRpt
#OUTPUT_DIR=`echo $INPUT | sed 's/\/S00.*.fastq.gz//g'`
#echo $OUTPUT_DIR
#make -f /gpfs/tgml/nin/exceRpt-master/exceRpt_smallRNA INPUT_FILE_PATH=/gpfs/tgml/reads/analysis_results/201709030_exceRpt/$INPUT OUTPUT_DIR=/gpfs/tgml/reads/analysis_results/201709030_exceRpt/${OUTPUT_DIR}_out MAIN_ORGANISM_GENOME_ID="hg38" DATABASE_PATH=/gpfs/tgml/nin/exceRpt/DATABASE ADAPTER_SEQ="none" MIN_READ_LENGTH=15 STAR_outFilterMatchNmin=15  ENDOGENOUS_LIB_PRIORITY=miRNA,tRNA,piRNA,gencode

#make -f /gpfs/tgml/nin/exceRpt-master/exceRpt_smallRNA INPUT_FILE_PATH=/gpfs/tgml/reads/analysis_results/201709030_exceRpt/$INPUT OUTPUT_DIR=/gpfs/tgml/reads/analysis_results/201709030_exceRpt/${OUTPUT_DIR}_out MAIN_ORGANISM_GENOME_ID="hg38" DATABASE_PATH=/gpfs/tgml/nin/exceRpt/DATABASE ADAPTER_SEQ="none" MIN_READ_LENGTH=15 STAR_outFilterMatchNmin=15  ENDOGENOUS_LIB_PRIORITY=miRNA,tRNA,piRNA,gencode
#out/tar/xvzf_exceRpt/exceRpt-4.3.2/exceRpt_smallRNA, out/tar/xvzf_exceRpt/exceRpt-4.3.2/mergePipelineRuns_functions.R, out/tar/xvzf_exceRpt/exceRpt-4.3.2/mergePipelineRuns.R

rule exceRpt:
    """
    Test: 
        out/exceRpt/smallRNA.gz/ln/updir/mw-el-cherif/inp/fastq/Run_245/S002433_batch-1_102_401-3_R1.fastq_CORE_RESULTS_v4.6.3.tgz
        out/exceRpt/smallRNA.gz/bedtools/bamtofastq/umi-tools/dedup/samtools/sam_to_bam_bai/bowtie2/se_ncrna-GRCh38-ensembl-r101/umi-tools/extract_se_NEXTflex-Small-RNA-seq-v3/cutadapt/se_-a_TGGAATTCTCGGGTGCCAAGG_--minimum-length_23/ln/updir/mw-el-cherif/inp/fastq/Run_245/S002350_batch-1_19_135-1_R1.fastq_CORE_RESULTS_v4.6.3.tgz
        out/exceRpt/smallRNA.gz/umi-tools/extract_se_NEXTflex-Small-RNA-seq-v3/cutadapt/se_-a_TGGAATTCTCGGGTGCCAAGG_--minimum-length_23/ln/updir/mw-el-cherif/inp/fastq/Run_245/S002350_batch-1_19_135-1_R1.fastq_CORE_RESULTS_v4.6.3.tgz
    """    
    input:
        fastq="out/{filler}.fastq{ext}",
        db="/gpfs/tgml/nin/exceRpt/DATABASE",
        exceRpt_smallRNA="/gpfs/tgml/nin/exceRpt-master/exceRpt_smallRNA",
        #exceRpt_smallRNA="out/unzip/exceRpt/exceRpt-master/exceRpt_smallRNA",
        #exceRpt_smallRNA="out/tar/xvzf_exceRpt/exceRpt-4.3.2/exceRpt_smallRNA"
    output:
        "out/{tool}{ext}/{filler}.fastq_CORE_RESULTS_v4.6.3.tgz"
        #"out/{tool}{ext}/{filler}.fastq_CORE_RESULTS_v5.0.0.tgz"
        #"out/{tool}{ext}/{filler}.fastq_CORE_RESULTS_v4.3.2.tgz"
    wildcard_constraints:
        ext="|.gz",
        tool="exceRpt/smallRNA"
    shell:
        """
        DIRNAME=`dirname {output}`
        make -f {input.exceRpt_smallRNA} INPUT_FILE_PATH={input.fastq} OUTPUT_DIR=$DIRNAME MAIN_ORGANISM_GENOME_ID="hg38" DATABASE_PATH={input.db} ADAPTER_SEQ="none" MIN_READ_LENGTH=15 STAR_outFilterMatchNmin=15 ENDOGENOUS_LIB_PRIORITY=miRNA,tRNA,piRNA,gencode
        """

rule exceRpt_merge_pipeline:
    """
    Test: 
        out/exceRpt/mergePipelineRuns/excerpt-test-S002350/done
    """
    input:
        excerpt_tgz = lambda wildcards: eval(config['ids'][wildcards.excerpt_tgz_list_id]),
        mergepipeline="/gpfs/tgml/nin/exceRpt-master/mergePipelineRuns.R",
    output:
        "out/exceRpt/mergePipelineRuns/{excerpt_tgz_list_id}/done"
    conda:
        "../../../../mw-lib/src/snakemake/envs/R_exceRpt.yaml"
    shell:
        """
        OUTDIR=out/exceRpt/mergePipelineRuns/{wildcards.excerpt_tgz_list_id}
        echo $OUTDIR
        mkdir -p $OUTDIR
        for file in {input.excerpt_tgz};
            do echo $file;
            #ln -srf $file $OUTDIR/;
            tar -zxvf $file -C $OUTDIR/;
        done
        echo $PWD
        cd $OUTDIR
        echo $PWD
        Rscript {input.mergepipeline} .
        touch done 
        """
