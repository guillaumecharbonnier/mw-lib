## TODO: Write pyatac_vplot_extra here and move rules below to mw-legacy

rule pyatac_vplot_lower_upper_flank_not_atac_plot_extra_bed:
    """
    Created:
        2017-06-07 15:07:53
    Aim:
        Produce vplots from sorted bam.
    Doc:
        https://nucleoatac.readthedocs.io/en/latest/pyatac/#vplot
        $ opt/miniconda/envs/nucleoatac/bin/pyatac vplot -h
        Command run:  opt/miniconda/envs/nucleoatac/bin/pyatac vplot -h
        pyatac version 0.3.1
        start run at: 2017-06-07 15:05
        usage: pyatac vplot [-h] --bed bed_file --bam bam_file [--out basename]
                            [--cores int] [--lower int] [--upper int] [--flank int]
                            [--scale] [--weight int] [--strand int] [--not_atac]
                            [--no_plot] [--plot_extra]
        
        optional arguments:
          -h, --help      show this help message and exit
        
        Required:
          Necessary arguments
        
          --bed bed_file  Positions around which to generate VPlot
          --bam bam_file  Accepts sorted BAM file
        
        General options:
        
          --out basename
          --cores int     Number of cores to use
        
        VMat option:
          Size, scaling of VPlot
        
          --lower int     lower limit on insert size
          --upper int     upper limit on insert size
          --flank int     how many bases on each side of site (or center of site) to
                          include
          --scale         Scale each site
          --weight int    column in which weight information is included
          --strand int    column in which strand information is included
          --not_atac      Don't use atac offsets
        
        Plot options:
          Make plots?
        
          --no_plot       Don't plot output
          --plot_extra    Make some extra plots
    
    Test:
        out/pyatac/vplot_lower-30_upper-200_flank-100_not_atac_plot_extra_bed-mm9-ss-s-wt-danpos-band2-maxfuz-55-minsmt-300/samtools/sort/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band1.Vplot.eps
    """
    input:
        bam="out/{filler}.bam",
        bai="out/{filler}.bam.bai",
        bed = lambda wildcards: eval(config['ids'][wildcards.bed_id])
    output:
        vplot=           "out/pyatac/vplot_lower-{lower}_upper-{upper}_flank-{flank}_not_atac_plot_extra_bed-{bed_id}/{filler}.Vplot.eps",
        insertsizes=     "out/pyatac/vplot_lower-{lower}_upper-{upper}_flank-{flank}_not_atac_plot_extra_bed-{bed_id}/{filler}.InsertSizes.eps",
        insertionprofile="out/pyatac/vplot_lower-{lower}_upper-{upper}_flank-{flank}_not_atac_plot_extra_bed-{bed_id}/{filler}.InsertionProfile.eps"
    params:
        outdir="out/pyatac/vplot_lower-{lower}_upper-{upper}_flank-{flank}_not_atac_plot_extra_bed-{bed_id}/{filler}"
    wildcard_constraints:
        lower="[0-9]+",
        upper="[0-9]+",
        flank="[0-9]+",
        bed_id="[a-zA-Z0-9-]+"
    conda:
        "../envs/nucleoatac.yaml"
    threads:
        16
    shell:
        """
        pyatac vplot \
            --bed {input.bed} \
            --bam {input.bam} \
            --out {params.outdir} \
            --cores {threads} \
            --lower {wildcards.lower} \
            --upper {wildcards.upper} \
            --flank {wildcards.flank} \
            --not_atac \
            --plot_extra
        """

rule pyatac_vplot_lower_upper_flank_scale_not_atac_plot_extra_bed:
    """
    Created:
        2017-06-08 14:07:51
    Aim:
        Produce vplots from sorted bam.
    Doc:
          --scale         Scale each site
    Test:
        out/pyatac/vplot_lower-30_upper-200_flank-100_not_atac_plot_extra_bed-mm9-ss-s-wt-danpos-band2-maxfuz-55-minsmt-300/samtools/sort/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm9/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band1.Vplot.eps
    """
    input:
        bam="out/{filler}.bam",
        bai="out/{filler}.bam.bai",
        bed = lambda wildcards: eval(config['ids'][wildcards.bed_id])
    output:
        vplot=           "out/pyatac/vplot_lower-{lower}_upper-{upper}_flank-{flank}_scale_not_atac_plot_extra_bed-{bed_id}/{filler}.Vplot.eps",
        insertsizes=     "out/pyatac/vplot_lower-{lower}_upper-{upper}_flank-{flank}_scale_not_atac_plot_extra_bed-{bed_id}/{filler}.InsertSizes.eps",
        insertionprofile="out/pyatac/vplot_lower-{lower}_upper-{upper}_flank-{flank}_scale_not_atac_plot_extra_bed-{bed_id}/{filler}.InsertionProfile.eps"
    params:
        outdir=          "out/pyatac/vplot_lower-{lower}_upper-{upper}_flank-{flank}_scale_not_atac_plot_extra_bed-{bed_id}/{filler}"
    wildcard_constraints:
        lower="[0-9]+",
        upper="[0-9]+",
        flank="[0-9]+",
        bed_id="[a-zA-Z0-9-]+"
    threads:
        16
    conda:
        "../envs/nucleoatac.yaml"
    shell:
        """
        pyatac vplot \
            --bed {input.bed} \
            --bam {input.bam} \
            --out {params.outdir} \
            --cores {threads} \
            --lower {wildcards.lower} \
            --upper {wildcards.upper} \
            --flank {wildcards.flank} \
            --scale \
            --not_atac \
            --plot_extra
        """

