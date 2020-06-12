rule dev_smoove:
    input:
        expand("out/imsindel/_fa-genome-GRCh38-r94-chr1/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_pe_fa-genome-GRCh38-r94-chr1/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/{sample}/1.out", sample=["Jurkat_SRR1509753_H3K27ac","Jurkat1_SRX2975115_H3K27ac","Jurkat2_SRR5789197_H3K27ac"])

rule smoove_single_call:
    """
    Created:
        2020-05-27 20:15:33
    Doc:
        https://github.com/brentp/smoove
    Warning:
        Should only work for paired-end reads
    Test:
        out/smoove/call_fa-genome-GRCh38/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_pe_fa-genome-GRCh38/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/1070_H3K27ac-smoove.genotyped.vcf.gz
    """
    input:
        bam="out/{filler}.bam",
        bai="out/{filler}.bam.bai",
        fa= lambda wildcards: eval(config['ids'][wildcards.fa_genome_id])
    output:
        vcf = "out/{tool}{extra}_{fa_genome_id}/{filler}-smoove.genotyped.vcf.gz"
    log:
              "out/{tool}{extra}_{fa_genome_id}/{filler}.log"
    benchmark:
              "out/{tool}{extra}_{fa_genome_id}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="smoove/call"
    threads: 2 #From doc: smoove can only parallelize up to 2 or 3 threads on a single-sample and it's most efficient to use 1 thread
    conda:
        "../envs/smoove.yaml"
    shell:
        """
        #smoove call -x --name my-cohort --exclude $bed --fasta $fasta -p $threads --genotype /path/to/*.bam
        #smoove call --outdir results-smoove/ --exclude $bed --name $sample --fasta $fasta -p {threads} --genotype {input.bam}
        SAMPLE=`basename {output.vcf} | sed 's/-smoove.genotyped.vcf.gz//'`
        OUTDIR=`dirname {output.vcf}`
        smoove call --outdir $OUTDIR --name $SAMPLE --fasta {input.fa} -p {threads} --genotype {input.bam} &> {log}
        """

