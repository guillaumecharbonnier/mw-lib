rule minimap2_ext:
    """
    Doc:
        https://lh3.github.io/minimap2/minimap2.html
        https://github.com/lh3/minimap2
    Test:
        out/minimap2/fa_fa-genome-GRCh38-ensembl-r100/megahit/pe/sickle/pe_-t_sanger_-q_30/ln/alias/sst/all_samples/fastq/PEER_Crispr_5-10_Chip_H3K27ac/final.contigs.sam
    Note:
        Defaut minimap setting do not map the contig containing the insertion to chr1. Testing other settings...

grep CACAGAAAGACGGTTAGGAAACGGTAACCCTACTTCCTGGCAGA out/meta/index_fa/gunzip/to-stdout/wget/ftp/ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa

grep TCTGCCAGGAAGTAGGGTTACCGTTTCCTAACCGTCTTTCTGTG out/meta/index_fa/gunzip/to-stdout/wget/ftp/ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa

Note: with blastn or dissimilar megablast I can see I have huge sections of the contig matching to nothing known, thus explaining why no match of the contif with bowtie2 or minimap

    """
    input:
        fa_to_align = "out/{filler}.{ext}",
        fa_genome = lambda wildcards: eval(mwconf['ids'][wildcards.fa_genome_id])
    output:
        sam = "out/{tool}{ext}{extra}_{fa_genome_id}/{filler}.sam"
    log:
              "out/{tool}{ext}{extra}_{fa_genome_id}/{filler}.log"
    benchmark:
              "out/{tool}{ext}{extra}_{fa_genome_id}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="minimap2/",
        ext="fa|fasta"
    threads: 2 #From doc: smoove can only parallelize up to 2 or 3 threads on a single-sample and it's most efficient to use 1 thread
    conda:
        "../envs/minimap2.yaml"
    shell:
        """
        minimap2 {params.extra} -a {input.fa_genome} {input.fa_to_align} > {output.sam} 2> {log}
        """

