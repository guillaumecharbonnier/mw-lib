
rule bedtools_makewindows_g_extra:
    """
    Created:
        2019-09-05 16:45:26
    Aim:
    Note:
        No doc online, look for command line
        Need -g (genome file) or -b (BED file) for interval source.
        Need -w (window size) or -n (number of windows).·
    Test:
        out/bedtools/makewindows_g_-w_100_-i_src/gunzip/to-stdout/rsync/ucsc/goldenPath/mm10/database/chromInfo.bed


        awk/convert_bed3_to_saf/bedtools/makewindows_g_-w_100/gunzip/to-stdout/rsync/ucsc/goldenPath/mm10/database/chromInfo.bed
    """
    input:
        chrominfo="out/{filler}.txt"
    output:
        bed="out/{tool}{extra}/{filler}.bed"
    log:
        "out/{tool}{extra}/{filler}.log"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="bedtools/makewindows_g"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools makewindows -g {input.chrominfo} {params.extra} > {output.bed} 2> {log}"
        

rule bedtools_makewindows_legacy:
    """
    Created:
        2017-06-19 14:42:19
    Aim:
    Test:
        out/bedtools/makewindows_w-200_i-src/ChromHMM/LearnModel_test5_numstates-20_assembly-mm10/spermatozoa_20_segments.bed
    """
    input:
        bed="out/{filler}.bed"
    output:
        bed="out/bedtools/makewindows_w-{w}_i-{i}/{filler}.bed"
    wildcard_constraints:
        w="[0-9]+",
        i="src|winnum|srcwinnum"
    conda:
        "../envs/bedtools.yaml"
    threads:
        1
    shell:
        """
        bedtools makewindows \
            -b {input.bed} \
            -w {wildcards.w} \
            -i {wildcards.i} > {output.bed}
        """


