rule ucsc_twoBitToFa:
    """
    Created:
        2024-10-24 10:43:51
    Aim:
        Convert 2bit to fa
    Test:
        out/ucsc/twoBitToFa/wget/https/hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa
    """
    input:
        twobit = "out/{filler}.2bit"
    output:
        fa = "out/ucsc/twoBitToFa/{filler}.fa"
    conda:
        "../envs/ucsc_twoBitToFa.yaml"
    shell:
        """
        twoBitToFa {input.twobit} {output.fa}
        """

