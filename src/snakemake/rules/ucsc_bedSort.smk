rule ucsc_bedSort:
    """
    Created:
        2018-03-06 18:31:08
    Aim:
	Sort bed-like files to be accepted by other UCSC tools like bedGraphToBigWig
    """
    input:
        bed = "out/{filler}.{ext}",
    output:
        bed = "out/ucsc/bedSort/{filler}.{ext}"
    conda:
        "../envs/ucsc_bedSort.yaml"
    shell:
        """
        bedSort {input.bed} {output.bed}
        """

