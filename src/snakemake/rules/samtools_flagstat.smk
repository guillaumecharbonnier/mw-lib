rule samtools_flagstat_in_input_directory:
    """
    Created:
        2017-04-11 14:45:26
    Aim:
        Give some stats about bam file. Used mainly for Capstarrseq analysis.
    """
    input:
        bam="out/{filler}.bam",
        bai="out/{filler}.bam.bai" # not sure the index is needed here.
    output:
        flagstat="out/{filler}.flagstat.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools flagstat {input.bam} > {output.flagstat}"
        
rule samtools_flagstat:
    """
    Created:
        2017-06-13 17:11:42
    Aim:
        Give some stats about bam file. Used mainly for Capstarrseq analysis.
    """
    input:
        bam="out/{filler}.bam",
        bai="out/{filler}.bam.bai" # not sure the index is needed here.
    output:
        flagstat="out/samtools/flagstat/{filler}.tsv"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools flagstat {input.bam} > {output.flagstat}"
