rule qualimap_bamqc:
    """
    Created:
        2019-02-15 15:36:47
    Doc:
        http://qualimap.bioinfo.cipf.es/doc_html/analysis.html#bam-qc
    Note:
        Bam have to be sorted first, else qualimap fails:
        java.lang.RuntimeException: The alignment file is unsorted.
    Test:
        out/qualimap/bamqc/samtools/sort/samtools/view_sam_to_bam/bowtie2/se_mm10/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR1202037/report.pdf
        out/qualimap/bamqc/samtools/sort/samtools/view_sam_to_bam/bowtie2/se_mm10/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR3126243/report.pdf
    """
    input:
        bam="out/{filler}.bam"
    output:
        pdf="out/qualimap/bamqc/{filler}/report.pdf"
    params:
        outdir="out/qualimap/bamqc/{filler}"
    conda:
        "../envs/qualimap.yaml"
    threads: 
        4
    shell:
        """
        qualimap bamqc -bam {input.bam} --paint-chromosome-limits -nt {threads} --genome-gc-distr mm9 -outformat PDF -outdir {params.outdir} -outfile report
        """
