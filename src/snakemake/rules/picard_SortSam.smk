rule picard_SortSam_sortOrder:
    """
    Created:
        2018-03-16 17:24:03
    Doc:
        https://gatk.broadinstitute.org/hc/en-us/articles/360036510732-SortSam-Picard
    Test:
        out/picard/SortSam_sortOrder-coordinate/samtools/view_bSh/bwa/samse_q-5_fa-GRCh38-Blueprint/gunzip/merge_lanes_nextseq500_se_raw/inp/fastq/run170/S001580_TH125_EC_K27ac-41046/TH125-EC-K27ac_S6.bam
    Note:
        queryname (Sorts according to the readname. This will place read-pairs and other derived
        reads (secondary and supplementary) adjacent to each other. Note that the readnames are
        compared lexicographically, even though they may include numbers. In paired reads, Read1
        sorts before Read2.)
        coordinate (Sorts primarily according to the SEQ and POS fields of the record. The
        sequence will sorted according to the order in the sequence dictionary, taken from from
        the header of the file. Within each reference sequence, the reads are sorted by the
        position. Unmapped reads whose mates are mapped will be placed near their mates. Unmapped
        read-pairs are placed after all the mapped reads and their mates.)
        duplicate (Sorts the reads so that duplicates reads are adjacent. Required that the
        mate-cigar (MC) tag is present. The resulting will be sorted by library, unclipped 5-prime
        position, orientation, and mate's unclipped 5-prime position.)

    """
    input:
        bam="out/{id}.bam",
    output:
        bam = "out/picard/SortSam_sortOrder-{sortOrder}/{id}.bam",
    wildcard_constraints:
        sortOrder="queryname|reads|coordinate|duplicate"
    conda:
        "../envs/picard.yaml"
    shell:
        """
        picard\
            SortSam\
            INPUT={input.bam}\
            OUTPUT={output.bam}\
            SORT_ORDER={wildcards.sortOrder}\
            VALIDATION_STRINGENCY=SILENT
    """

