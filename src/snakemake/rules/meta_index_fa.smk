rule meta_index_fa:
    """
    Created:
        2019-12-19 16:04:30
    Aim:
        produce both a dict and a fai index file for a fasta reference file.
        Dict file is required by some gatk functions.
    Note:
        Check samtoools can index a gz file.
        picard can handle it.
    Test:
        out/meta/index_genome_fa/cat/assembly_ensembl/GRCh38.fa
    """
    input:
        fasta="out/{filler}.{ext}"
        #fasta="out/{filler}.fasta"
    output:
        fasta="out/meta/index_{ext}/{filler}.{ext}",
        faidx="out/meta/index_{ext}/{filler}.{ext}.fai",
        dict="out/meta/index_{ext}/{filler}.dict"
    threads: 1
    wildcard_constraints:
        ext="fa|fasta|fa.gz|fasta.gz"
    conda:
        "../envs/picard.yaml"
    shell:
        """
        ln -srf {input.fasta} {output.fasta}
        picard CreateSequenceDictionary R={output.fasta} O={output.dict}
        samtools faidx {output.fasta}
        """
