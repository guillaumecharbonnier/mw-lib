rule picard_CreateSequenceDictionary:
    """
    Created:
        2019-12-19 15:11:36
    Test:
        out/picard/CreateSequenceDictionary_fa/cat/assembly_ensembl/GRCh38.fa
        out/picard/CreateSequenceDictionary_fa/bwa/aggregate_fa_genome_and_sam_contigs_with_inserts_fa-genome-hg19-main-chr/bowtie2/se_fasta_--rfg_1,1_-k_1_-q_-f_hg19/edena/assembling_-d_20_-c_20_-minCoverage_5/edena/overlapping/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_H3K27ac_contigs.dict
    """
    input:
        fasta="out/{filler}.{ext}"
        #fasta="out/{filler}.fasta"
    output:
        fasta="out/picard/CreateSequenceDictionary_{ext}/{filler}.{ext}",
        dict="out/picard/CreateSequenceDictionary_{ext}/{filler}.dict"
    threads: 1
    wildcard_constraints:
        ext="fa|fasta|fa.gz|fasta.gz"
    conda:
        "../envs/picard.yaml"
    shell:
        """
        ln -srf {input.fasta} {output.fasta}
        picard CreateSequenceDictionary R={input.fasta} O={output.dict}
        """

