rule bedtools_makewindows:
    """
    Created:
        2017-06-19 14:42:19
    Aim:
    Note:
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


