rule samtools_faidx:
    """
    Created:
        2017-05-16 11:37:45
    Aim:
         Index reference sequence in the FASTA format or extract subsequence from indexed reference sequence. If no region is specified, faidx will index the file and create <ref.fasta>.fai on the disk. If regions are specified, the subsequences will be retrieved and printed to stdout in the FASTA format.

         The input file can be compressed in the BGZF format.

         The sequences in the input file should all have different names. If they do not, indexing will emit a warning about duplicate sequences and retrieval will only produce subsequences from the first sequence with the duplicated name. 
    Test:
        out/samtools/faidx/gunzip/to-stdout/rsync/ucsc/goldenPath/hg38/chromosomes/chr19.fa
    """
    input:
        fasta = "out/{filler}.{ext}"
    output:
        fasta = "out/samtools/faidx/{filler}.{ext}",
        fai   = "out/samtools/faidx/{filler}.{ext}.fai"
    wildcard_constraints:
        ext="fasta|fa"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        ln -srf {input.fasta} {output.fasta}
        samtools faidx {output.fasta}
        """
