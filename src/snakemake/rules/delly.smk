rule dev_delly:
    input:
        expand("out/imsindel/_fa-genome-GRCh38-r94-chr1/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_pe_fa-genome-GRCh38-r94-chr1/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/{sample}/1.out", sample=["Jurkat_SRR1509753_H3K27ac","Jurkat1_SRX2975115_H3K27ac","Jurkat2_SRR5789197_H3K27ac"])

rule delly_call:
    """
    Created:
        2020-05-27 19:34:19
    Doc:
        https://github.com/dellytools/delly
    Warning:
        Should only work for paired-end reads
        Should also work on single-end reads:
        https://github.com/dellytools/delly/issues/183#issuecomment-591450099

    Test:
        out/delly/call_fa-genome-GRCh38/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_pe_fa-genome-GRCh38/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/1070_H3K27ac.bcf
        out/delly/call_-i_fa-genome-GRCh38/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_pe_fa-genome-GRCh38/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/1070_H3K27ac.bcf

    """
    input:
        bam="out/{filler}.bam",
        bai="out/{filler}.bam.bai",
        fa= lambda wildcards: eval(config['ids'][wildcards.fa_genome_id])
    output:
        bcf = "out/{tool}{extra}_{fa_genome_id}/{filler}.bcf"
    log:
              "out/{tool}{extra}_{fa_genome_id}/{filler}.log"
    benchmark:
              "out/{tool}{extra}_{fa_genome_id}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="delly/call"
    threads: MAX_THREADS
    conda:
        "../envs/delly.yaml"
    shell:
        """
        export OMP_NUM_THREADS={threads}
        delly call -g {input.fa} -o {output.bcf} {input.bam} &> {log}
        """

