rule genometools_bed_to_gff3:
    """
    Created:
        2018-03-12 16:21:03
    Aim:
        Convert bed files to gff3, mainly to use with ChIPSeqSpike
    Test:
        out/genometools/bed_to_gff3/ln/alias/microarray-upreg-ko-nut-threshold-1.5.gff
        out/gt/bed_to_gff3/awk/chromInfo_to_bed3/awk/extract_main_chr/gunzip/to-stdout/rsync/ucsc/goldenPath/hg38/database/chromInfo.gff
    """
    input:
        bed="out/{filler}.bed",
    output:
        gff="out/gt/bed_to_gff3/{filler}.gff"
    conda:
        "../envs/gt.yaml"
    shell:
        "gt bed_to_gff3 -o {output.gff} {input.bed}"
