rule bedtools_subtract_prepare_exclude_for_shuffle:
    """
    """
    input:
        non_unmapable="out/awk/get_non_unmapable/ucsc/bigWigToBedGraph/rsync/ucsc/gbdb/mm9/bbi/crgMapabilityAlign75mer.bedGraph",
        chromInfo="out/mysql/ucsc/mm9/chromInfo_main.tsv"
    output:
        unmapable="out/bedtools/subtract/prepare_exclude_for_shuffle/ucsc/bigWigToBedGraph/rsync/ucsc/gbdb/mm9/bbi/crgMapabilityAlign75mer.bed"
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        cat {input.chromInfo} | awk 'BEGIN {{OFS="\\t"}}; {{print $1, 0, $2}}' | \
        bedtools subtract -a - -b {input.non_unmapable} > {output.unmapable}
        """

