rule bedtools_merge:
    """
    Created:
        2017-03-13 11:16:53
    Aim:
        Merge overlapping feature in a bed file. First use is to merge WT and KO peaks for H4K5ac in order to do see if the ratio is centered on 0 for controls
    Note:
        bedtools merge requires that you presort your data by chromosome and then by start position (e.g., sort -k1,1 -k2,2n in.bed > in.sorted.bed for BED files).
    Suggested preceding rule:
        sort_coordinates_bed
    Test:
        out/bedtools/merge/sort/coordinates-bed/cat/mm10_nut_wt_ko_h4k5ac_peaks.bed
        out/bedtools/merge/sort/coordinates-bed/cat/mm9-P5424-K27ac-WT-KO-IK-peaks.bed
        out/bedtools/merge/sort/coordinates-bed/cat/hg38-fastkd1-run229-peaks.bed
        out/bedtools/merge/sort/_-k1,1_-k2,2n/bedtools/bamtobed/samtools/sort/samtools/view_sam_to_bam/bowtie2/se_mm10/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR1202037.bed
        out/bedtools/merge/sort/coordinates-bed/bedtools/bamtobed/samtools/sort/samtools/view_sam_to_bam/bowtie2/se_mm10/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR1202037.bed
    """
    input:
        bed="out/{filler}.bed"
    output:
        bed="out/{tool}{extra}/{filler}.bed"
    wildcard_constraints:
        tool="bedtools/merge"
    params:
        extra = params_extra
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools merge {params.extra} -i {input.bed} > {output.bed}"


