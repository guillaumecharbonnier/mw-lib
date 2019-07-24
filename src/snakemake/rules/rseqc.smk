rule rseqc_geneBody_coverage:
    """
    Created:
        2016-11-3 11h04
    Aim:
    Doc:
        http://rseqc.sourceforge.net/#genebody-coverage-py
    Test:
        
    """
    input:
        bam = "out/{filler}.bam",
        bai = "out/{filler}.bam.bai",
        bed = "annotation/processed/feature/{index}/rseqc/housekeeping.bed"
    output:
        pdf       = "out/rseqc/geneBody_coverage/{filler}.geneBodyCoverage.curves.pdf"
    params:
        outprefix = "out/rseqc/geneBody_coverage/{filler}"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "geneBody_coverage.py -i {input.bam} -r {input.bed} -o {params.outprefix}"

def input_bam_rseqc_geneBody_coverage_bam_list(wildcards):
    # Getting variables from Snakemake wildcards.
    bam_list_key=wildcards['bam_list_key']
    index=wildcards['index']
    
    if bam_list_key == 'nut':
        prefix="out/star/pe_" + index + "/sickle/pe_-t_sanger_-q_20/fastq/nut_rnaseq/"
        paths=expand("{prefix}{sample}.bam", prefix=prefix, sample=["Nut-P-WT","Nut-P-KO","Nut-R-WT","Nut-R-KO"])

    elif bam_list_key == 'h2al2':
        paths_p_r_c=expand("out/star/pe_GRCm38_outFilterMultimapNmax-1000/sickle/pe_-t_sanger_-q_20/merge_lanes/run176/RNA-{stage}-H2AL2-{condition}-Rep{replicate}.bam",   
            stage=["P","R","C"],  
            condition=["WT","KO"],   
            replicate=["1","2","3"])

        paths_s_r1=expand("out/star/pe_GRCm38_outFilterMultimapNmax-1000/sickle/pe_-t_sanger_-q_20/merge_lanes/run176/RNA-S-H2AL2-{condition}-Rep1.bam",
            condition=["WT","KO"])

        paths_s_r2=expand("out/star/pe_GRCm38_outFilterMultimapNmax-1000/sickle/pe_-t_sanger_-q_20/merge_lanes/run176/RNA-S-H2AL2-WT-Rep2.bam")
       
        paths = paths_p_r_c + paths_s_r1 + paths_s_r2

    else:
        print('Error: Add paths in function input_bam_rseqc_geneBody_coverage_bam_list')

    return(paths)


def input_bai_rseqc_geneBody_coverage_bam_list(wildcards):
    # Getting variables from Snakemake wildcards.
    bam_list_key=wildcards['bam_list_key']
    index=wildcards['index']
    
    if bam_list_key == 'nut':
        prefix="out/star/pe_" + index + "/sickle/pe_-t_sanger_-q_20/fastq/nut_rnaseq/"
        paths=expand("{prefix}{sample}.bam.bai", prefix=prefix, sample=["Nut-P-WT","Nut-P-KO","Nut-R-WT","Nut-R-KO"])

    elif bam_list_key == 'h2al2':
        paths_p_r_c=expand("out/star/pe_GRCm38_outFilterMultimapNmax-1000/sickle/pe_-t_sanger_-q_20/merge_lanes/run176/RNA-{stage}-H2AL2-{condition}-Rep{replicate}.bam.bai",   
            stage=["P","R","C"],  
            condition=["WT","KO"],   
            replicate=["1","2","3"])

        paths_s_r1=expand("out/star/pe_GRCm38_outFilterMultimapNmax-1000/sickle/pe_-t_sanger_-q_20/merge_lanes/run176/RNA-S-H2AL2-{condition}-Rep1.bam.bai",
            condition=["WT","KO"])

        paths_s_r2=expand("out/star/pe_GRCm38_outFilterMultimapNmax-1000/sickle/pe_-t_sanger_-q_20/merge_lanes/run176/RNA-S-H2AL2-WT-Rep2.bam.bai")
       
        paths = paths_p_r_c + paths_s_r1 + paths_s_r2

    else:
        print('Error: Add paths in function input_bam_rseqc_geneBody_coverage_bam_list')

    return(paths)


rule rseqc_geneBody_coverage_bam_list:
    """
    Created:
        2016-11-3 14h29
    Modified:
        2017-03-28 15:23:41 - Updated paths for H2AL2 analysis.
    Note:
        This should not work properly outside H2AL2 case for the moment after rework.
    Test:
        out/rseqc/geneBody_coverage_bam_list/mm10/h2al2.geneBodyCoverage.curves.pdf
    """
    input:
        bam=input_bam_rseqc_geneBody_coverage_bam_list,
        bai=input_bai_rseqc_geneBody_coverage_bam_list,
        #bed="annotation/processed/feature/{index}/rseqc/housekeeping.bed"
        bed="out/sed/remove_chr/awk/extract_main_chr/gunzip/wget/sourceforge_rseqc/BED/Mouse_Mus_musculus/mm10.HouseKeepingGenes.bed"
    output:
        pdf="out/rseqc/geneBody_coverage_bam_list/{index}/{bam_list_key}.geneBodyCoverage.curves.pdf"
    params:
        outprefix="out/rseqc/geneBody_coverage_bam_list/{index}/{bam_list_key}"
    conda:
        "../envs/rseqc.yaml"
    shell:
        """
        BAM_LIST=`echo "{input.bam}" | tr ' ' ','`
        geneBody_coverage.py -i $BAM_LIST -r {input.bed} -o {params.outprefix}
        """

rule rseqc_infer_experiment:
    """
    Created:
        2019-02-05 23:25:41
    Doc:
        http://rseqc.sourceforge.net/#infer-experiment-py
    Test:
        out/rseqc/infer_experiment/bed-hg38-GENCODE-knownGene/ln/updir/mw/inp/bam/hg19_RNA-Seq_thymus/ISP.bam.txt
        Does not work for this bam and I do not know why:
        out/rseqc/infer_experiment/bed-hg38-GENCODE-knownGene/out/star/pe_staridx-GRCh38_gtf-GRCh38-ensembl/sickle/pe_-t_sanger_-q_30/sra-tools/fastq-dump_pe/SRR2040277.bam.txt
    """
    input:
        sam_or_bam="out/{filler}.{sam_or_bam}",
        bed = lambda wildcards: config['ids'][wildcards.bed_id]
    output:
        txt = "out/rseqc/infer_experiment/{bed_id}/{filler}.{sam_or_bam}.txt"
    log:
              "out/rseqc/infer_experiment/{bed_id}/{filler}.{sam_or_bam}.log"
    benchmark:
              "out/rseqc/infer_experiment/{bed_id}/{filler}.{sam_or_bam}.benchmark.tsv"
    wildcard_constraints:
        sam_or_bam='sam|bam',
        bed_id='bed-[a-zA-Z0-9-]+'
    conda:
        "../envs/rseqc.yaml"
    shell:
        "infer_experiment.py -r {input.bed} -i {input.sam_or_bam} > {output.txt} 2> {log}"

    

