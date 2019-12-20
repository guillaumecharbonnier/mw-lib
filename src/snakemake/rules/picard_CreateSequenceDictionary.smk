rule picard_CreateSequenceDictionary:
    """
    Created:
        2019-12-19 15:11:36
    Test:
        out/picard/CreateSequenceDictionary_fa/cat/assembly_ensembl/GRCh38.fa
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

