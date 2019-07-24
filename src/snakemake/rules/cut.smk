rule cut_extra:
    """
    Created:
        2017-10-11 10:47:54
    Aim:
        Cut based on column field
    Test:
        out/cut/_-f1-3/awk/select_subpopulations_from_bed_lmin-4000/ln/broadPeak_to_bed/gunzip/to-stdout/wget/ftp/hgdownload.cse.ucsc.edu/apache/htdocs/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneA549H3k04me3Dex100nmPk.bed
    """
    input:
        txt="out/{filler}"
    output:
        txt="out/{tool}{extra}/{filler}"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="cut/"
    conda:
        "../envs/coreutils.yaml"
    shell:
        "cut {params.extra} {input.txt} > {output.txt}"
