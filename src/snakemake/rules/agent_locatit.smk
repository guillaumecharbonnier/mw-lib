rule agent_locatit_mbc:
    """
    Aim:
        This version of Locatit requires MBC file, and is for example required for RNA-seq in XT-HS2 mode (Single-strand consensus)
    Doc:
        https://www.agilent.com/cs/library/software/public/AGeNT%20ReadMe.pdf
    Note:
        Aligned reads should be sorted by read names and not by coordinates.
    Test:
        out/samtools/index/samtools/sort/agent/locatit_mbc_-i_-R/picard/SortSam_sortOrder-queryname/star/pe_fastq.gz_to_bam_staridx-GRCh38-ensembl_gtf-GRCh38-ensembl/agent/trim_-v2/ln/updir/mw/inp/fastq/2021_RNAseq_NECKER_spicuglia/fastq/MOLT4_S60.bam
    """
    input:
        bam = "out/{filler_align}{filler_trim}.bam",
        mbc = "out/{filler_trim}.txt.gz",
        agent = "out/agent/agent/agent.sh",
        locatit = "out/agent/agent/lib/locatit-2.0.5.jar"
    output:
        bam = "out/{tool}{extra}/{filler_align}{filler_trim}.bam"
    log:
        "out/{tool}{extra}/{filler_align}{filler_trim}.log"
    benchmark:
        "out/{tool}{extra}/{filler_align}{filler_trim}.benchmark.tsv"
    params:
        extra = params_extra,
        memory="120G" # arbitrary value, need refinement, 64 is enough from most samples, but largest samples need more.
    threads:
        MAX_THREADS
    wildcard_constraints:
        tool = "agent/locatit_mbc",
        filler_trim = "agent/trim.*"
    conda:
        "../envs/agent.yaml"
    shell:
        """
        java -Xmx{params.memory} -jar {input.locatit} {params.extra} -o {output.bam} {input.bam} {input.mbc} &> {log}
        """
