"""
Created:
    2017-05-10 16:01:41
Aim:
    Here are gathered rules with used multiple chained core tools.
"""

rule sh_extract_gfold_gene_id_to_bed3_lower_than_threshold:
    """
    Created:
        2017-03-21 14:38:28
    Aim:
        Convert gfold gene_id that contain coordinates to true bed.
    Test:
        "out/sh/extract_gfold_gene_id_to_bed3_lower_than_threshold-2/r/sort_gfold_diff_by_gfold_then_log2ratio/gfold/diff/gfold/count_bed-h4k5ac-peaks/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/merge_lanes/run140_run141/H4K5ac-Nut-WT_VS_H4K5ac-Nut-KO_diff_ascend.bed" 
    """
    input:
        tsv="out/{filler}.tsv"
    output:
        bed="out/sh/extract_gfold_gene_id_to_bed3_lower_than_threshold{threshold}/{filler}.bed"
    wildcard_constraints:
        threshold="[0-9-]+"
    shell:
        """
        awk '$3 < {wildcards.threshold}' {input.tsv} | cut -f1 | tr -d '"' | tr '-' '\t' > {output.bed}
        """

rule sh_danpos_xls_to_bed3:
    """
    Modified:
        2017-05-10 15:47:04
    Aim:
        Remove header and unwanted columns to make bed3.
    """   
    input:
        xls="out/{filler}.xls"
    output:
        bed="out/sh/danpos_xls_to_bed3/{filler}.bed"
    shell: 
        """
        tail -n+2 {input.xls} | cut -f1-3 > {output.bed}
        """


