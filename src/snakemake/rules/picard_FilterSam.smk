
rule picard_FilterSamReads_extra:
    """
    Test:
        out/picard/FilterSamReads_JoseDavidContamination/ln/alias/sst/all_samples/GRCh38/bam/ATOJ7_CB548.bam
    """
    input:
        bam="out/{filler}.bam",
        #bai="out/{filler}.bam.bai" # bai may be not needed for this function. Testing this assumption on 2020-02-11 17:27:49
    output:
        bam =     "out/{tool}{extra}/{filler}.bam"
    params:
        tmp_readnames = "out/{tool}{extra}/{filler}.readnames.txt",
        extra = params_extra
    wildcard_constraints:
        tool = "picard/FilterSamReads"
    conda:
        "../envs/picard.yaml"
    shell:
        """
        samtools view {input.bam} | awk '($3 == "chr14" && $9 == 121 && $4 == 22448604)' | cut -f1 > {params.tmp_readnames}

        # java -jar picard.jar 
        picard FilterSamReads \
            I={input.bam} \
            O={output.bam} \
            READ_LIST_FILE={params.tmp_readnames} \
            FILTER=excludeReadList

        """
