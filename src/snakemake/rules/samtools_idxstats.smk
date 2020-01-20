rule samtools_idxstats:
    """
    Created:
        2017-05-16 11:30:54
    Aim:
        Retrieve and print stats in the index file corresponding to the input file. Before calling idxstats, the input BAM file must be indexed by samtools index.
        The output is TAB-delimited with each line consisting of reference sequence name, sequence length, # mapped reads and # unmapped reads. It is written to stdout.
    Note:
        'idxstat' needed by multiqc : https://multiqc.info/docs/#idxstats
    Test:
        out/samtools/idxstats/samtools/index/samtools/sort/samtools/view_sam_to_bam_-q_30/bowtie2/se_mm10/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR3126243.idxstat.tsv
        
        out/samtools/idxstats/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem2_se_bwa-index-hg19-main-chr-and-contigs-with-inserts-from-T11C-H3K27ac/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_H3K27ac.idxstat.tsv
    """
    input:
        bam="out/{filler}.bam",
        bai="out/{filler}.bam.bai"
    output:
        tsv="out/samtools/idxstats/{filler}.idxstat.tsv"
    conda:
        "../envs/samtools.yaml"
    threads:
        1
    shell:
        "samtools idxstats {input.bam} > {output.tsv}"
