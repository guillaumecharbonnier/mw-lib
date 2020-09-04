rule dev_imsindel:
    input:
        expand("out/imsindel/_fa-genome-GRCh38-r94-chr1/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_pe_fa-genome-GRCh38-r94-chr1/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/{sample}/1.out", sample=["Jurkat_SRR1509753_H3K27ac","Jurkat1_SRX2975115_H3K27ac","Jurkat2_SRR5789197_H3K27ac"])

rule dev_mindthegap_megahit:
    input:
        mindthegap_jurkat_se=expand("out/MindTheGap/fill/MindTheGap/find_se_fa-genome-GRCh38-ensembl-r100/gunzip/to-stdout/{sickle}ln/alias/sst/all_samples/fastq/{sample}.insertions.fasta", sickle=["","sickle/se_-t_sanger_-q_30/"], sample=['Jurkat1_SRR1603650_H3K27ac','Jurkat_SRR1057274_H3K27ac']),
        mindthegap_jurkat_pe=expand("out/MindTheGap/fill/MindTheGap/find_pe_fa-genome-GRCh38-ensembl-r100/gunzip/to-stdout/{sickle}ln/alias/sst/all_samples/fastq/{sample}.insertions.fasta", sickle=["","sickle/pe_-t_sanger_-q_30/"], sample=["Jurkat_SRR1509753_H3K27ac","Jurkat1_SRX2975115_H3K27ac","Jurkat2_SRR5789197_H3K27ac"]),
        mindthegap_tall_se=expand("out/MindTheGap/fill/MindTheGap/find_se_fa-genome-GRCh38-ensembl-r100/gunzip/to-stdout/{sickle}ln/alias/sst/all_samples/fastq/{sample}.insertions.fasta", sickle=["","sickle/se_-t_sanger_-q_30/"], sample=[x.strip() for x in open("../mw-tall/src/snakemake/lists/tall_chip_samples_se.txt","r")]),
        mindthegap_tall_pe=expand("out/MindTheGap/fill/MindTheGap/find_pe_fa-genome-GRCh38-ensembl-r100/gunzip/to-stdout/{sickle}ln/alias/sst/all_samples/fastq/{sample}.insertions.fasta", sickle=["","sickle/pe_-t_sanger_-q_30/"], sample=[x.strip() for x in open("../mw-tall/src/snakemake/lists/tall_chip_samples_pe.txt","r")]),
        megahit_jurkat_se=expand("out/samtools/idxstats/samtools/index/picard/MarkDuplicates_REMOVE_DUPLICATES=true_VALIDATION_STRINGENCY=SILENT/samtools/sort/samtools/view_sam_to_bam_-q_30/bowtie2/build_flanked_indel_megahit_and_align_se_GRCh38-ensembl-r100/{sickle}ln/alias/sst/all_samples/fastq/{sample}.idxstat.tsv", sickle=["","sickle/se_-t_sanger_-q_30/"], sample=['Jurkat1_SRR1603650_H3K27ac','Jurkat_SRR1057274_H3K27ac']),
        megahit_jurkat_pe=expand("out/samtools/idxstats/samtools/index/picard/MarkDuplicates_REMOVE_DUPLICATES=true_VALIDATION_STRINGENCY=SILENT/samtools/sort/samtools/view_sam_to_bam_-q_30/bowtie2/build_flanked_indel_megahit_and_align_pe_GRCh38-ensembl-r100/{sickle}ln/alias/sst/all_samples/fastq/{sample}.idxstat.tsv", sickle=["","sickle/se_-t_sanger_-q_30/"], sample=["Jurkat_SRR1509753_H3K27ac","Jurkat1_SRX2975115_H3K27ac","Jurkat2_SRR5789197_H3K27ac"]),
        megahit_tall_se=expand("out/samtools/idxstats/samtools/index/picard/MarkDuplicates_REMOVE_DUPLICATES=true_VALIDATION_STRINGENCY=SILENT/samtools/sort/samtools/view_sam_to_bam_-q_30/bowtie2/build_flanked_indel_megahit_and_align_se_GRCh38-ensembl-r100/{sickle}ln/alias/sst/all_samples/fastq/{sample}.idxstat.tsv", sickle=["","sickle/se_-t_sanger_-q_30/"], sample=[x.strip() for x in open("../mw-tall/src/snakemake/lists/tall_chip_samples_se.txt","r")]),
        megahit_tall_pe=expand("out/samtools/idxstats/samtools/index/picard/MarkDuplicates_REMOVE_DUPLICATES=true_VALIDATION_STRINGENCY=SILENT/samtools/sort/samtools/view_sam_to_bam_-q_30/bowtie2/build_flanked_indel_megahit_and_align_pe_GRCh38-ensembl-r100/{sickle}ln/alias/sst/all_samples/fastq/{sample}.idxstat.tsv", sickle=["","sickle/se_-t_sanger_-q_30/"], sample=[x.strip() for x in open("../mw-tall/src/snakemake/lists/tall_chip_samples_pe.txt","r")])

rule dev_filterpeaks:
    input:
        expand("out/bedtools/intersect_-sorted_-v_-b_bed-GRCh38-H3K27ac-thymocytes-peaks/sort/coordinates-bed/ln/alias/sst/all_samples/GRCh38/bed/broad/{sample}_peaks.bed", sample = [x.strip() for x in open("../mw-tall/src/snakemake/lists/tall_chip_samples_se.txt","r")])

rule imsindel:
    """
    Test:
        out/imsindel/_fa-genome-GRCh38-r94-chr1/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_se_fa-genome-GRCh38-r94-chr1/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/Jurkat_SRR1057274_H3K27ac.vcf
    """
    input:
        #imsindel="../IMSindel/bin/imsindel",
        bam="out/{filler}.bam",
        fa= lambda wildcards: eval(config['ids'][wildcards.fa_genome_id])
    output:
        vcf="out/{tool}{extra}_{fa_genome_id}/{filler}/1.out"
    log:
        "out/{tool}{extra}_{fa_genome_id}/{filler}/log"
    benchmark:
        "out/{tool}{extra}_{fa_genome_id}/{filler}/benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="imsindel/"
    threads:
        MAX_THREADS
    conda:
        "../envs/imsindel.yaml"
    shell:
        """
        imsindel --thread {threads} --bam {input.bam} --chr 1 --indelsize 10000 --outd `dirname {output}` --reffa {input.fa} &> {log}
        """

