rule deepTools_plotEnrichment:
    """
    Created:
        2017-10-05 10:39:03
    Note:
        plotEnrichment seems to take the whole path as label and regionLabel leading to error because of too long string in plot. Thus labels are given based on file name.
    Test:
        out/deepTools/plotEnrichment_bed-mm10-test-srr-peaks_bam-mm10-test-srr.eps
    """
    input:
        bam_list = lambda wildcards: eval(config['ids'][wildcards.bam_list_id]),
        bai_list = lambda wildcards: [path + '.bai' for path in eval(config['ids'][wildcards.bam_list_id])],
        bed_list = lambda wildcards: eval(config['ids'][wildcards.bed_list_id])
    output:
        png_pdf_svg_eps_plotly = "out/{tool}{extra}_{bed_list_id}_{bam_list_id}.{png_pdf_svg_eps_plotly}"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "deepTools/plotEnrichment",
        png_pdf_svg_eps_plotly = "png|pdf|svg|eps|plotly"
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        LABELS=""
        for i in {input.bam_list}
        do
            LABELS="$LABELS `basename $i | sed 's/.bam//g'`"
        done
        echo $LABELS

        REGION_LABELS=""
        for i in {input.bed_list}
        do
            REGION_LABELS="$REGION_LABELS `basename $i | sed 's/.bed//g'`"
        done
        echo $REGION_LABELS

        plotEnrichment {params}\
            --bamfiles {input.bam_list} \
            --BED {input.bed_list} \
            --plotFile {output.png_pdf_svg_eps_plotly} \
            --labels $LABELS \
            --regionLabels $REGION_LABELS
        """

