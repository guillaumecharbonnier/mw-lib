rule pasha_input:
    """
    out/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/ln/alias/fastq/run215/H4K5ac-Nut-WT.bam out/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/ln/alias/fastq/run215/H4K5ac-Nut-WT.bam out/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_dm6/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/ln/alias/fastq/run215/H4K5ac-Nut-WT.bam out/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/ln/alias/fastq/run215/Input-Nut-WT.bam out/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/ln/alias/fastq/run215/H4K5ac-Nut-KO.bam out/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_dm6/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/ln/alias/fastq/run215/H4K5ac-Nut-KO.bam out/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/ln/alias/fastq/run215/Input-Nut-KO.bam
    
    
    out/ucsc/wigToBigWig_mm10_clip/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/ln/alias/fastq/run215/pasha/result_PE/AllReads/WIGfs_H4K5ac-Nut-WT_mergedReads_elPairsAndEst153_AThr3_bin50.bw
    ,out/ucsc/wigToBigWig_mm10_clip/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/ln/alias/fastq/run215/pasha/result_PE/AllReads/WIGfs_Input-Nut-WT_mergedReads_elPairsAndEst152_AThr4_bin50.bw
    
    
    
    out/ucsc/wigToBigWig_mm10_clip/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/ln/alias/fastq/run215/pasha/result_PE/AllReads/WIGfs_H4K5ac-Nut-KO_mergedReads_elPairsAndEst152_AThr3_bin50.bw,out/ucsc/wigToBigWig_mm10_clip/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/ln/alias/fastq/run215/pasha/result_PE/AllReads/WIGfs_Input-Nut-KO_mergedReads_elPairsAndEst153_AThr3_bin50.bw
    """

rule pasha_pe:
    """
    Created:
        2018-08-19 19:26:36
    Aim:
        Produce wig to be converted to bigwig to use with ChIPSeqSpike because the ones produced by deepTools bamCoverage produces bugs.
    Test:
        out/pasha/ln/alias/experiments/ChIP-Seq_Spike-in_REH/GRCh38/run219-Input-gev.bam
    """
    input:
        #Rscript="opt/miniconda/envs/pasha/bin/Rscript",
        script="src/r/script/pasha.R",
        bam="out/{filler}/{basename}.bam"
    output:
        bam="out/pasha/l-{se_pe}/{filler}/{basename}.bam",
        wig="out/pasha/l-{se_pe}/{filler}/{basename}.wig"
    benchmark:
        "log/snakemake/benchmark/pasha/l-{se_pe}/{filler}/{basename}.log"
    conda:
        "../envs/r_pasha.yaml"
    threads:
        1
    shell:
        """
        ln -f {input.bam} {output.bam}
        Rscript {input.script} -b {output.bam} -t {threads} -w {output.wig} -l {wildcards.se_pe}
        
        OUTWIG=`find out/pasha/l-{wildcards.se_pe} -wholename '*/{wildcards.filler}/*/AllReads/WIGfs_{wildcards.basename}.bam_*.wig'`
        echo $OUTWIG | wc -w 
        # Should always be 1 and print warning if not.
        ln -f $OUTWIG {output.wig}
        """

