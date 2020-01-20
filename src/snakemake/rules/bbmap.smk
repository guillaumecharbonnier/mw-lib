rule bbmap_se:
    """
    Created:
        2019-12-02 01:30:30
    Aim:
        First used to select subpopulations from bam.
        Trying this:
        https://www.biostars.org/p/225198/
    Note:
    Test:
        out/bbmap/se_fastq.gz_fa-genome-hg19-contigs-with-insert-in-T11C-H3K27ac/ln/alias/sst/all_samples/fastq/T11C_H3K27ac.sam
        Old:
            out/bbmap/se_fa-genome-GRCh38-r94-chr1/ln/alias/sst/all_samples/fastq/Jurkat_SRR1057274_H3K27ac.sam.gz
    """
    input:
        fq="out/{filler}.{ext}",
        fa= lambda wildcards: eval(config['ids'][wildcards.fa_genome_id])
    output:
        sam="out/{tool}_{ext}{extra}_{fa_genome_id}/{filler}.sam"
    log:
        "out/{tool}_{ext}{extra}_{fa_genome_id}/{filler}.log"
    benchmark:
        "out/{tool}_{ext}{extra}_{fa_genome_id}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    threads:
        MAX_THREADS
    wildcard_constraints:
        tool = "bbmap/se",
        ext = "fastq.gz|fasta" #ADD MORE HERE IF NEEDED
    conda:
        "../envs/bbmap.yaml"
    shell:
        "bbmap.sh ref={input.fa} in={input.fq} out={output.sam} {params.extra} t={threads} maxindel=100k &> {log}"
