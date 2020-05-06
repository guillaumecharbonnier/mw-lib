rule ucsc_bedToBigBed:
    """
    Created:
        2019-09-06 18:38:32
    Note:
        Need a sorted bed:
        sort -k1,1 -k2,2n unsorted.bed > input.bed
        Need to remove lines starting with "track"

        bed9
    Test:
        out/ucsc/bedToBigBed_chrominfo-hg19/ChromHMM/MakeBrowserFiles/ChromHMM/MakeSegments_our_merged_thymic_samples_into_chromdet_original_paper/CD34_11_segments.bb

    """
    input:
        bed = "out/{filler}.bed",
        chromInfo = lambda wildcards: eval(config['ids'][wildcards.chrominfo_id])
    output:
        bb="out/{tool}{extra}_{chrominfo_id}/{filler}.bb"
    log:
        "out/{tool}{extra}_{chrominfo_id}/{filler}.log"
    benchmark:
        "out/{tool}{extra}_{chrominfo_id}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="ucsc/bedToBigBed"
    conda:
        "../envs/ucsc.yaml"
    shell:
        "bedToBigBed {params.extra} {input.bed} {input.chromInfo} {output.bb} &> {log}"
