rule texlive_pdfcrop:
    """
    Created:
        2016-03-04 11:07:00
    Modified:
        2017-06-08 16:16:24 - Move from bedtools_jaccard to its own file.
    Aim:
        It is quite hard to define the size of the pdf output containing a table created in R with gridExtra package. This rule is a dirty hack that allow to hardcode very large width and height in "format_jaccard_danpos_peaks.R" and then use pdfcrop to remove all the blank not used for each table.
        Also allow to crop Vplots produced by pyatac.
    Test:
        out/texlive/pdfcrop/texlive/epstopdf/pyatac/vplot_lower-30_upper-200_flank-300_not_atac_plot_extra_bed-mm9-ss-s-wt-danpos-band2-maxfuz-55-minsmt-300/samtools/sort/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band1.Vplot.pdf
    """
    input:
        pdf="out/{filler}.pdf",
    output:
        pdf="out/texlive/pdfcrop/{filler}.pdf"
    conda:
        "../envs/texlive_selected.yaml"
    threads:
        1
    shell:
        """
        pdfcrop {input.pdf} {output.pdf}
        """

rule texlive_pdfjam:
    """
    Created:
        2018-11-05 19:08:00
    Modified:
        2017-06-08 16:16:24 - Move from bedtools_jaccard to its own file.
    Aim:
    Note:
        Does not currently work because:
            pdfjam ERROR: LaTeX package pdfpages.sty is not installed
        I will just write beamer slides for my current needs.
    Test
        out/texlive/pdfjam_--nup_2x1_--landscape/test.pdf
    """
    input:
        pdf = lambda wildcards: eval(mwconf['ids'][wildcards.pdf_list_id])
    output:
        pdf = "out/{tool}{extra}/{pdf_list_id}"
    params:
        extra = params_extra
    conda:
        "../envs/texlive_selected.yaml"
    wildcard_constraints:
        tool = "texlive/pdfjam"
    threads:
        1
    shell:
        """
        pdfjam {params.extra} --outfile {output.pdf} {input.pdf}
        """
