rule jvarkit_bam2wig:
    """
    Doc:
        http://lindenb.github.io/jvarkit/Bam2Wig.html
    Test:
        
    """
    input:
        "out/{filler}.bam"
    output:
        "out/jvarkit/bam2wig/{filler}.wig"
    conda:
        "../envs/jvarkit_bam2wig.yaml"
    shell:
        "bam2wig.sh -o {output} {input}"
