rule picard_FastqToSam_pe_dev_for_locatit_without_mbc:
    """
    Aim:
        Added as a component for Agilent Agent pipeline, but actually unneeded, hence uncomplete and untested.
    """
    input:
        fq1 = "out/{filler}_1.fastq.gz",
        fq2 = "out/{filler}_2.fastq.gz"
    output:
        bam = "out/{tool}{extra}/{filler}.bam"
    wildcard_constraints:
        tool = "picard/FastqToSam"
    shell:
        """
        java -jar picard.jar FastqToSam \
            F1={input.fq1} \
            F2={input.fq2} \
            O={output.bam} \
            SM=sample001 \
            RG=rg0013
        """