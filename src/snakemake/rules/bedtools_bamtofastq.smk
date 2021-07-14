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
        tool="bedtools/bamtofastq_se"
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

rule bedtools_bamtofastq_pe_extra:
    """
    Created:
    Aim:
    Note:
        Input should be sorted by read name
        samtools sort -n aln.bam aln.qsort
    Test:
        out/bedtools/bamtofastq_pe/agent/locatit_mbc_-i_-R/picard/SortSam_sortOrder-queryname/star/pe_fastq.gz_to_bam_staridx-GRCh38-ensembl_gtf-GRCh38-ensembl/agent/trim_-v2/ln/updir/mw/inp/fastq/2021_RNAseq_NECKER_spicuglia/fastq/MOLT4_S60_1.fastq.gz
    """
    input:
        bam = "out/{filler}.bam"
    output:
        fq1 = "out/{tool}{extra}/{filler}_1.fastq.gz",
        fq2 = "out/{tool}{extra}/{filler}_2.fastq.gz"
    log:
        "out/{tool}{extra}/{filler}.log"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="bedtools/bamtofastq_pe"
    conda:
        "../envs/bedtools.yaml"
    threads:
        1
    shell:
        """
        FQ1=`echo {output.fq1} | sed 's/.gz$//'`
        FQ2=`echo {output.fq2} | sed 's/.gz$//'`
        bedtools bamtofastq {params.extra} -i {input.bam} -fq $FQ1 -fq2 $FQ2 2> {log}
        gzip $FQ1 2>> {log}
        gzip $FQ2 2>> {log}
        """

