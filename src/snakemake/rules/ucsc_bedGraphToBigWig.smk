rule ucsc_bedGraphToBigWig:
    """
    Created:
        2018-03-06 18:31:08
    Note:
        The input bedGraph file must be sorted, use the unix sort command:
          sort -k1,1 -k2,2n unsorted.bedGraph > sorted.bedGraph
          options:
            -blockSize=N - Number of items to bundle in r-tree.  Default 256
            -itemsPerSlot=N - Number of data points bundled at lowest level. Default 1024
            -unc - If set, do not use compression.
    Test:
        out/ucsc/bedGraphToBigWig_mm10/sort/coordinates_bed/bedtools/genomecov_bga_ibam_g-mm10/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/ln/alias/fastq/run207/H4K5ac-Nut-WT.bw
    """
    input:
        BedGraph = "out/{filler}.BedGraph",
        chromInfo = lambda wildcards: config['ids'][wildcards.chrominfo_id]
    output:
        bw="out/ucsc/bedGraphToBigWig_{chrominfo_id}/{filler}.bw"
    conda:
        "../envs/ucsc.yaml"
    shell:
        "bedGraphToBigWig {input.BedGraph} {input.chromInfo} {output.bw}"
