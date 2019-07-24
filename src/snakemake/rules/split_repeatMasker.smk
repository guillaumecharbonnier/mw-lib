rule split_repeatMasker:
    """
    Created: 2016-03-31 17h04
    This rule split the repeatMasker bed file into individual bed file for each repeat type.
    
    We select the 4th column in repeatMasker file which contain the name of the repeat type, then we keep only one of each type for the loop. Anchors tabulations are mandatory for grep because a lot of repeat type names are variations of the same pattern, e.g. searching L1ME3 will also retrieve regions called L1ME3A, L1ME3B and L1ME3C.
    
    Usage:
    expand("annotation/processed/feature/{index}/split_repeatMasker/", index=["mm9","mm10"])
    """
    input:
        bed="annotation/input/feature/{index}/repeatMasker.bed"
    output:
        dir="annotation/processed/feature/{index}/split_repeatMasker/"
    threads: 1
    shell:"""
    # Defining a tab variable because \t does not work well in grep. 
    TAB=`echo -e "\t"`
    
    for REPEAT_TYPE in `cut -f4 {input.bed} | sort | uniq`
    do
        REPEAT_NO_SPECIAL_CHAR=`echo ${{REPEAT_TYPE}} | tr '[()]' '_'`
        grep "${{TAB}}${{REPEAT_TYPE}}${{TAB}}" {input.bed} > {output.dir}${{REPEAT_NO_SPECIAL_CHAR}}.bed
    done
    """

rule awk_extract_feature_in_repeatMasker:
    """
    Created: 2017-01-09 11h56 - Adapted to new pattern concept and added wanted feature in output.

    This rule split the repeatMasker bed file into individual bed file for each repeat type.
    
    We select the 4th column in repeatMasker file which contain the name of the repeat type, then we keep only one of each type for the loop. Anchors tabulations are mandatory for grep because a lot of repeat type names are variations of the same pattern, e.g. searching L1ME3 will also retrieve regions called L1ME3A, L1ME3B and L1ME3C.
    
    Test:
        "out/awk/extract_feature_in_repeatMasker/mm9/AT_rich.bed"
    """
    input:
        bed="inp/annotation/input/feature/{index}/repeatMasker.bed"
    output:
        bed="out/awk/extract_feature_in_repeatMasker/{index}/{feature}.bed"
    threads:
        1
    wildcard_constraints:
        feature="[a-zA-Z0-9_-]+"
    shell:
        """
        awk '$4=="{wildcards.feature}"' {input.bed} > {output.bed}
        """

