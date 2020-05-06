rule danpos_dtriple_extra:
    """
    Created:
        2019-01-22 17:11:46
    Doc:
          -m , --paired         set to 1 if the input data is mate-pair (paired-end) reads.  (default: 0)
    Test:
        out/danpos/dtriple/ln/alias/experiments/hg38_H3K27ac_LAL-T/CD34.wig -fn
        bed-mm9-ss-c-wtA
        bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run113_run119_run124/s1-MNS-R-WT_s2-MNS-R-KO.bed

        out/awk/extract_main_chr/awk/convert_bedpe_to_bed6_insert_size/bedtools/bamtobed_bedpe/samtools/sort_n/samtools/merge_two_samples/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run113_run119_run124/s1-MNS-R-WT_s2-MNS-R-KO.bed

out/samtools/merge_three_runs/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/ln/updir/mw-sk/inp/fastq/tgml/run113_run119_run124/MNS-R-WT.bam


out/samtools/merge_three_runs/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/pe_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run113_run119_run124/MNS-R-WT.bam

../mw-sk/inp/fastq/tgml/run119/

    """
    input:
        bam="out/{filler}.bam"#,
        #python="opt/miniconda/envs/danpos/bin/python",
        #danpos="opt/danpos-2.2.2/danpos.py"
    output:
        wig =         "out/{tool}{extra}/{filler}/danpos.smooth.wig",
        ln_wig =      "out/{tool}{extra}/{filler}.wig",
        xls_positions="out/{tool}{extra}/{filler}/danpos.smooth.positions.xls",
        xls_peaks =   "out/{tool}{extra}/{filler}/danpos.smooth.peaks.xls",
        xls_regions = "out/{tool}{extra}/{filler}/danpos.smooth.regions.xls"
    log:
                      "out/{tool}{extra}/{filler}.log"
    benchmark:
                      "out/{tool}{extra}/{filler}.benchmark.tsv"
    params:
        extra = params_extra,
        ln_bam =      "out/{tool}{extra}/{filler}/danpos.bam",
        outdir =      "out/{tool}{extra}/{filler}/"
    wildcard_constraints:
        tool="danpos/dtriple"
    threads:
        1
    priority:
        2
    conda:
        "../envs/danpos.yaml"
    shell:
        """
        # Danpos include input file name in the output directory, which is ugly.
        # Link is used here so Danpos see clean input file.
        # Note: hardlink is bad with snakemake because when the hardlink is done on the bam there is a 'touch' on file, thus every rules using the bam have to be done again.
        WDIR=`pwd`
        
        (ln -s --force $WDIR/{input.bam} $WDIR/{params.ln_bam}
        cd {params.outdir}

        INPUT_SAMPLE=`basename {params.ln_bam}`

        danpos.py dtriple {params.extra} --out . $INPUT_SAMPLE

        # Danpos put output files in './result/pooled/' when given '.' as '--out'.
        mv result/pooled/* ./
        rm -rf result
        ln $WDIR/{output.wig} $WDIR/{output.ln_wig} ) &> {log}
        """

