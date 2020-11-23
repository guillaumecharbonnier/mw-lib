rule bedtools_bamtofastq_extra:
    """
    Created:
    Aim:
    Note:
    Test:
        out/bedtools/bamtofastq/umi-tools/dedup/samtools/sam_to_bam_bai/bowtie2/se_ncrna-GRCh38-ensembl-r101/umi-tools/extract_se_NEXTflex-Small-RNA-seq-v3/cutadapt/se_-a_TGGAATTCTCGGGTGCCAAGG_--minimum-length_23/ln/updir/mw-el-cherif/inp/fastq/Run_245/S002350_batch-1_19_135-1_R1.fastq.gz

    """
    input:
        bam="out/{filler}.bam"
    output:
        "out/{tool}{extra}/{filler}.fastq.gz"
    log:
        "out/{tool}{extra}/{filler}.log"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="bedtools/bamtofastq"
    conda:
        "../envs/bedtools.yaml"
    threads:
        1
    shell:
        """
        TMP=`echo {output} | sed 's/.gz$//'`
        bedtools bamtofastq {params.extra} -i {input.bam} -fq $TMP 2> {log}
        gzip $TMP 2>> {log}
        """
