rule bedtools_complement:
    """
    Created:
        2017-04-12 14:56:27
    Aim:
        bedtools complement returns all intervals in a genome that are not covered by at least one interval in the input BED/GFF/VCF file. 
        Used for salva to get a bed file for deepTools blacklist.
    Test:
        out/bedtools/complement_hg38/sort/coordinates_bed/awk/extract_main_chr/crossmap/bed_hg19_to_hg38/input/bed/salva/Regions_capture_SE_Starr-seq_Alex_HG19_Merged_cellLine_specific_SE.bed
    """
    input:
        chrominfo="out/sort/coordinates_bed/awk/extract_main_chr/awk/extract_two_first_columns/gunzip/rsync/ucsc/goldenPath/{index}/database/chromInfo.txt",
        bed="out/{filler}.bed",
    output:
        bed="out/bedtools/complement_{index}/{filler}.bed"
    log:
            "out/bedtools/complement_{index}/{filler}.log"
    benchmark:
            "out/bedtools/complement_{index}/{filler}.benchmark.tsv"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools complement -i {input.bed} -g {input.chrominfo} > {output.bed} 2> {log}"
