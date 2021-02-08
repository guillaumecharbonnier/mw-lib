rule bedtools_closest_extra:
    """
    Created:
        2017-03-15 17:52:33
    Modified:
        2018-09-19 15:06:24 - First time I try this 'extra' pattern. May be good to try if it could work for other "duplicated" rules.
    Aim:
        -d  In addition to the closest feature in B, report its distance to A as an extra column. The reported distance for overlapping features will be 0.
    Doc:
        https://bedtools.readthedocs.io/en/latest/content/tools/closest.html
    Note:
        bedtools closest requires that all input files are presorted data by chromosome and then by start position (e.g., sort -k1,1 -k2,2n in.bed > in.sorted.bed for BED files).
    Test:
        out/bedtools/closest_-D_ref_-id_-k_1_-t_first_bed-hg38-ensembl-r93-tss-protein-coding-TR-IG/sort/bed/ln/alias/experiments/hg38_atac_peaks_from_wilson/CD34-EC.bed
        out/bedtools/closest_-D_ref_-iu_-k_1_-t_first_bed-hg38-ensembl-r93-tss-protein-coding-TR-IG/sort/bed/ln/alias/experiments/hg38_atac_peaks_from_wilson/CD34-EC.bed
    """
    input:
        a = "out/{filler}.bed",
        b = lambda wildcards: eval(mwconf['ids'][wildcards.bed_id])
    output:
        bed = "out/{tool}{extra}_{bed_id}/{filler}.bed"
    log:
              "out/{tool}{extra}_{bed_id}/{filler}.log"
    benchmark:
              "out/{tool}{extra}_{bed_id}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool  = "bedtools/closest",
        bed_id = "bed-[a-zA-Z0-9-]+"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools closest -a {input.a} -b {input.b} {params.extra} > {output.bed} 2> {log}"
