rule ucsc_bigWigToBedGraph:
    """
    Created:
        2016-12-07 16h20 - First attemp to convert bw from crgmappability.
    
    Example:
        "out/ucsc/bigWigToBedGraph/rsync/ucsc/gbdb/mm9/bbi/crgMapabilityAlign75mer.bedGraph"
    """
    input:
        bw="out/{filler}.bw"
    output:
        bedGraph="out/ucsc/bigWigToBedGraph/{filler}.bedGraph"
    conda:
        "../envs/ucsc.yaml"
    shell:
        "bigWigToBedGraph {input.bw} {output.bedGraph}"
