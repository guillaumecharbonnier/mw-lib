rule paste_great_input_and_output:
    """
    Created:
        2018-03-13 16:59:16
    Aim:
        Paste bed files given as input for GREAT and the output with distance to nearby genes.
    Test:
        out/paste/great_input_and_output_hg19/crossmap/bed_hg38_to_hg19/bedtools/intersect_thymus_stage_common_peaks_hg38_v3/not_common_subset_3.bed

    Note:
    sed -e 's/[^0-9]/ /g' -e 's/^ *//g' -e 's/ *$//g' out/sed/remove_first_and_last_line/doc/great/

    """
    input:
        bed="out/{filler}.bed",
        txt="out/coreutils/extract_min_dist_from_great/sed/remove_first_and_last_line/doc/great/{assembly}/{filler}.txt"
    output:
        bed="out/paste/great_input_and_output_{assembly}/{filler}.bed"
    shell:
        """
        paste {input.bed} {input.txt} > {output.bed}
        """
