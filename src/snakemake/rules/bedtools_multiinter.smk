rule bedtools_multiinter_bed:
    """
    Created:
        2018-04-23 15:31:26
    Aim:
        Default bedtools multiinter rule
    Note:
        Requires that each interval file is sorted by chrom/start.
    Test:
        out/bedtools/multiinter/bed-h3k27ac-supe1000-minus-blacklist-minus-1250-TSS-1250.bed
        out/bedtools/multiinter_-header/bed-hg38-sequencing-summary-h3k27ac-2019-10-26.bed
    """
    input:
        bed_list = lambda wildcards: eval(mwconf['ids'][wildcards.bed_list_id])
    output:
        bed_multiinter = "out/{tool}{extra}/{bed_list_id}.bed"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="bedtools/multiinter"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools multiinter {params.extra} -i {input.bed_list} > {output.bed_multiinter}"

