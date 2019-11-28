rule ucsc_faToTwoBit:
    """
    Created:
        2019-11-26 22:21:43
    Test:
        out/ucsc/faToTwoBit/gunzip/to-stdout/wget/http/hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr22.2bit
    """
    input:
        fa = "out/{filler}.fa",
    output:
        twobit="out/ucsc/faToTwoBit/{filler}.2bit"
    conda:
        "../envs/ucsc.yaml"
    shell:
        "faToTwoBit {input.fa} {output.twobit}"
