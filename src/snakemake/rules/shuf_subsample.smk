rule shuf_subsample:
    """
    Created:
        2016-12-20 20h01
    Aim:
        Created to do random sampling in a bed file but could be applied to any text file. Tested as an alternative to sort. Do benchmark and see which one is the fastest.
    Note:
        wc -l out/awk/filter_bedpe_by_size/min40_max80/samtools/view/bampe_to_bed3_nodup/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126/PSK-SC-KO.bed
        69771672 out/awk/filter_bedpe_by_size/min40_max80/samtools/view/bampe_to_bed3_nodup/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126/PSK-SC-KO.bed

        The number of line in this file should be used to subsample PSK-SC-WT

    Test:
        "out/shuf/subsample69771672/awk/filter_bedpe_by_size/min40_max80/samtools/view/bampe_to_bed3_nodup/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126/PSK-SC-WT.bed"

    Error for large sampling:
        shuf -n 125549155 out/awk/filter_bedpe_by_size/min40_max80/samtools/view/bampe_to_bed3_nodup/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126/PSK-SC-WT.bed > out/shuf/subsample125549155/awk/filter_bedpe_by_size/min40_max80/samtools/view/bampe_to_bed3_nodup/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126/PSK-SC-WT.bed

        shuf: memory exhausted
    
    Test:
        out/shuf/subsample1000/awk/extract_feature_in_repeatMasker/mm9/AT_rich.bed
    """
    input:
        txt="out/{filler}"
    output:
        txt="out/shuf/subsample{n}/{filler}"
    wildcard_constraints:
        n="[0-9]+"
    shell:"""
    INPUT_N_LINES=`wc -l {input.txt}`
    echo $INPUT_N_LINES

    shuf -n {wildcards.n} {input.txt} > {output.txt}
    """
