rule gtftk_bed_to_gtf_featureType_source:
    """
    Created:
        2017-02-21 16:43:17
    Warning:
        Probably a removed function as it is not available in
        https://pygtftk.readthedocs.io/en/latest/
    Aim:
        Convert a bed to gtf. Useful for example to  use featureCounts on unconventionnal features.
    Test:
        out/gtftk/bed_to_gtf_featureTypetest_sourcegc_workflow/awk/extract_xls_coordinates_to_bed6/macs2/callpeak_broad_BAMPE_mm_SPMRTrue/bam/mm9/merge_run140_run141/H4K5ac-Nut-WT_over_Input-Nut-WT_peaks.gtf
    """
    input:
        bed="out/{filler}.bed"
    output:
        gtf="out/gtftk/bed_to_gtf_featureType{featureType}_source{source}/{filler}.gtf"
    wildcard_constraints:
        featureType="transript|test",
        source="Unknown|gc_workflow"
    conda:
        "../envs/pygtftk.yaml"
    shell:
        """
        gtftk bed_to_gtf \
            -i {input.bed} \
            -o {output.gtf} \
            -t {wildcards.featureType} \
            -s {wildcards.source}
        """


