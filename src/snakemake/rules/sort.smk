rule sort_extra:
    """
    Created:
        2019-01-30 12:35:22
    Test:
        out/sort/bed/crossmap/hg38_to_hg19/ln/alias/experiments/hg38_atac_peaks_from_wilson/CD34-EC.bed
        out/sort/_-k1,1_-k2,2n/crossmap/hg38_to_hg19/ln/alias/experiments/hg38_atac_peaks_from_wilson/CD34-EC.bed
    """
    input:
        "out/{filler}"
    output:
        "out/{tool}{extra}/{filler}"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="sort/"
    shell:
        "LC_ALL=C; sort {params.extra} {input} > {output}"

ruleorder: sort_coordinates_bedlike_legacy > sort_extra
ruleorder: sort_unique_coordinates_bed_legacy > sort_extra
ruleorder: sort_coordinates_rseg_bed_legacy > sort_extra

rule sort_unique_coordinates_bed_legacy:
    """
    Created:
        2017-01-24 16h36
    Aim:
        Rule to remove redundant coordinates (chr/start/end) for features. Very annoying when you have multiple time the same TSS in heatmaps.
    Test:
        out/sort/unique_coordinates_bed/input/annotation/feature/debugged_ranked_florent/mm9_5000_features_upreg.bed
    """
    input:
        "out/{id}.bed"
    output:
        "out/sort/unique_coordinates_bed/{id}.bed"
    shell:
        """
        sort -u -k1,3 {input} > {output}
        """

rule sort_coordinates_bedlike_legacy:
    """
    Created:
        2017-03-13 14:31:09
    Modified:
        2017-04-12 15:29:14 - Removed .bed suffix because I want to be able to use it also to sort chrominfo.
        2019-01-30 15:35:28 - Changed to legacy, you can use 'out/sort/bed' now instead
    Aim:
        Sort data by chromosome and then by start position. Required before to use bedtools merge. Or to use memory-efficient algorithm in bedtools coverage.
    Test:
    """
    input:
        "out/{filler}"
    output:
        "out/sort/coordinates_bed/{filler}"
    shell:
        """
        sort -k1,1 -k2,2n {input} > {output}
        """

rule sort_coordinates_rseg_bed_legacy:
    """
    Created:
        2017-11-02 13:25:11
    Aim:
        Sort data as requested as input by RSEG Test:
    """
    input:
        "out/{filler}"
    output:
        "out/sort/coordinates_rseg_bed/{filler}"
    shell:
        """
        export LC_ALL=C
        sort -k1,1 -k2,2n -k3,3n -k6,6 -o {output} {input}
        """

