rule dqRNASeq:
    """
    Created:
        2020-02-10 17:39:38
    Aim:
        Testing dqRNASeq to compare with UMI-tools dedup
    Doc:
        https://github.com/e-hutchins/dqRNASeq
    """
    input:
        script="../dqRNASeq/dqRNASeq.sh",
        stl="inp/STL96.txt",
        bam="out/star/pe_fastq.gz_to_bam_standard_staridx-GRCh38-ensembl_gtf-GRCh38-ensembl/{filler}.bam",
        fq1="out/{filler}_1.fastq.gz",
        fq2="out/{filler}_2.fastq.gz"
    output:
    shell:
        """
        {input.script} -b {input.bam} -f {input.fq1} -r {input.fq2} -t {input.stl}
        """
