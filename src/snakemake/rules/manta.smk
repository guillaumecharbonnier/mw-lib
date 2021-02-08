rule dev_manta:
    input:
        expand("out/imsindel/_fa-genome-GRCh38-r94-chr1/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_pe_fa-genome-GRCh38-r94-chr1/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/{sample}/1.out", sample=["Jurkat_SRR1509753_H3K27ac","Jurkat1_SRX2975115_H3K27ac","Jurkat2_SRR5789197_H3K27ac"])

rule manta_single_call:
    """
    Created:
        2020-05-27 20:57:11
    Doc:
        https://github.com/Illumina/manta
    Warning:
        Should only work for paired-end reads
    Test:
        out/manta/call_fa-genome-GRCh38/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_pe_fa-genome-GRCh38/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/1070_H3K27ac.log
    """
    input:
        bam="out/{filler}.bam",
        bai="out/{filler}.bam.bai",
        fa= lambda wildcards: eval(mwconf['ids'][wildcards.fa_genome_id])
    #output:
    #    done = "out/{tool}{extra}_{fa_genome_id}/{filler}.done"
    # ADD this later:
    # https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md#structural-variant-predictions
    log:
              "out/{tool}{extra}_{fa_genome_id}/{filler}.log"
    benchmark:
              "out/{tool}{extra}_{fa_genome_id}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="manta/call"
    threads: 2 #From doc: smoove can only parallelize up to 2 or 3 threads on a single-sample and it's most efficient to use 1 thread
    conda:
        "../envs/manta.yaml"
    shell:
        """
        configManta.py --bam {input.bam} --referenceFasta {input.fa} --runDir `dirname {log}`
        """

