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
        out/ucsc/bedGraphToBigWig_chrominfo-hg19/igvtools/tdftobedgraph/ln/updir/mw-tall/inp/GSE60104/GSM1464990_20140211_733.spikein.hg19.bedgraph.bw
    """
    input:
        bedGraph = "out/{filler}.bedGraph",
        chromInfo = lambda wildcards: eval(mwconf['ids'][wildcards.chrominfo_id])
    output:
        bw="out/ucsc/bedGraphToBigWig_{chrominfo_id}/{filler}.bw"
    conda:
        "../envs/ucsc.yaml"
    shell:
        """
        grep -v "track type=" {input.bedGraph} | \
            sort -k1,1 -k2,2n > {output}.tmp.bedGraph
            bedGraphToBigWig {output}.tmp.bedGraph {input.chromInfo} {output.bw}
        rm {output}.tmp.bedGraph
        """

rule ucsc_BedGraphToBigWig:
    """
    Created:
        2023-02-06 18:31:08
    Aim:
        Added to keep track of legacy code where bedGraph was written BedGraph in input
	This rule should be moved in mw-legacy when I have time
    Note:
        The input bedGraph file must be sorted, use the unix sort command:
          sort -k1,1 -k2,2n unsorted.bedGraph > sorted.bedGraph
          options:
            -blockSize=N - Number of items to bundle in r-tree.  Default 256
            -itemsPerSlot=N - Number of data points bundled at lowest level. Default 1024
            -unc - If set, do not use compression.
    Test:
        out/ucsc/bedGraphToBigWig_mm10/sort/coordinates_bed/bedtools/genomecov_bga_ibam_g-mm10/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/ln/alias/fastq/run207/H4K5ac-Nut-WT.bw
        out/ucsc/bedGraphToBigWig_chrominfo-hg19/igvtools/tdftobedgraph/ln/updir/mw-tall/inp/GSE60104/GSM1464990_20140211_733.spikein.hg19.bedgraph.bw
    """
    input:
        BedGraph = "out/{filler}.BedGraph",
        chromInfo = lambda wildcards: eval(mwconf['ids'][wildcards.chrominfo_id])
    output:
        bw="out/ucsc/BedGraphToBigWig_{chrominfo_id}/{filler}.bw"
    conda:
        "../envs/ucsc.yaml"
    shell:
        """
        bedGraphToBigWig {input.BedGraph} {input.chromInfo} {output.bw}
        """
