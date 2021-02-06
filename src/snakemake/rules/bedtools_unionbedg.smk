rule bedtools_unionbedg:
    """
    Created:
        2019-10-24 22:38:31
    Aim:
        bedtools unionbedg combines multiple BEDGRAPH files into a single file
        such that one can directly compare coverage (and other text-values
        such as genotypes) across multiple sample
    Doc:
        https://bedtools.readthedocs.io/en/latest/content/tools/unionbedg.html
    Note:
        NAMES is required because the header argument does not work properly
        to label each input bed file.
    Test:
        out/bedtools/unionbedg/bed-hg19-segments-merged-thymic-samples.bed
        out/bedtools/unionbedg_-header/bed-hg19-segments-blueprint-tall-samples.bg
        out/bedtools/unionbedg_-header/bed-hg19-5-states-collapsed-segments-thymopoiesis-tall-samples.bg
    """
    input:
        bed_list = lambda wildcards: eval(mwconf['ids'][wildcards.bg_list_id])
    output:
        bed = "out/{tool}{extra}/{bg_list_id}.bg"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="bedtools/unionbedg"
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        NAMES=`basename -a {input} | sed 's/\.[a-z]*$//'`
        echo $NAMES
        bedtools unionbedg {params.extra} -i {input.bed_list} -names $NAMES> {output.bed}
        """

