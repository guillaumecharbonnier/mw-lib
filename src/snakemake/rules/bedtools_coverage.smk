rule bedtools_coverage:
    """
    Created:
        2016-08-24 14h59
    Aim:
        The bedtools coverage tool computes both the depth and breadth of coverage of features in file B on the features in file A. For example, bedtools coverage can compute the coverage of sequence alignments (file B) across 1 kilobase (arbitrary) windows (file A) tiling a genome of interest. One advantage that bedtools coverage offers is that it not only counts the number of features that overlap an interval in file A, it also computes the fraction of bases in the interval in A that were overlapped by one or more features. Thus, bedtools coverage also computes the breadth of coverage observed for each interval in A.
    Doc:
        https://bedtools.readthedocs.io/en/latest/content/tools/coverage.html
        Option	Description
        -a	BAM/BED/GFF/VCF file “A”. Each feature in A is compared to B in search of overlaps. Use “stdin” if passing A with a UNIX pipe.
        -b	One or more BAM/BED/GFF/VCF file(s) “B”. Use “stdin” if passing B with a UNIX pipe. NEW!!!: -b may be followed with multiple databases and/or wildcard (*) character(s).
        -abam	BAM file A. Each BAM alignment in A is compared to B in search of overlaps. Use “stdin” if passing A with a UNIX pipe: For example: samtools view -b <BAM> | bedtools intersect -abam stdin -b genes.bed. Note: no longer necessary after version 2.19.0
        -hist	
        Report a histogram of coverage for each feature in A as well as a summary histogram for _all_ features in A.
        Output (tab delimited) after each feature in A:
        1) depth
        2) # bases at depth
        3) size of A
        4) % of A at depth
        -d	Report the depth at each position in each A feature. Positions reported are one based. Each position and depth follow the complete A feature.
        -counts	Only report the count of overlaps, don’t compute fraction, etc. Restricted by -f and -r.
        -f	Minimum overlap required as a fraction of A. Default is 1E-9 (i.e. 1bp).
        -F	Minimum overlap required as a fraction of B. Default is 1E-9 (i.e., 1bp).
        -r	Require that the fraction of overlap be reciprocal for A and B. In other words, if -f is 0.90 and -r is used, this requires that B overlap at least 90% of A and that A also overlaps at least 90% of B.
        -e	Require that the minimum fraction be satisfied for A _OR_ B. In other words, if -e is used with -f 0.90 and -F 0.10 this requires that either 90% of A is covered OR 10% of B is covered. Without -e, both fractions would have to be satisfied.
        -s	Force “strandedness”. That is, only report hits in B that overlap A on the same strand. By default, overlaps are reported without respect to strand.
        -S	Require different strandedness. That is, only report hits in B that overlap A on the _opposite_ strand. By default, overlaps are reported without respect to strand.
        -split	Treat “split” BAM (i.e., having an “N” CIGAR operation) or BED12 entries as distinct BED intervals.
        -sorted	For very large B files, invoke a “sweeping” algorithm that requires position-sorted (e.g., sort -k1,1 -k2,2n for BED files) input. When using -sorted, memory usage remains low even for very large files.
        -g	Specify a genome file the defines the expected chromosome order in the input files for use with the -sorted option.
        -header	Print the header from the A file prior to results.
        -sortout	When using multiple databases (-b), sort the output DB hits for each record.
        -nobuf	Disable buffered output. Using this option will cause each line of output to be printed as it is generated, rather than saved in a buffer. This will make printing large output files noticeably slower, but can be useful in conjunction with other software tools and scripts that need to process one line of bedtools output at a time.
        -iobuf	Follow with desired integer size of read buffer. Optional suffixes K/M/G supported. Note: currently has no effect with compressed files.
    Note:
        Draft, maybe deprecated for my current need in comparison with deepTools multbigwigsummary + plotCorrelation
    Test:
        

    """
    input:
        a="out/{filler}.{ext_a}",
        b=lambda wildcards: config['ids'][wildcards.b_id]
    output:
        "out/{tool}{extra}_{b_id}/{filler}.bed"
    log:
        "out/{tool}{extra}_{b_id}/{filler}.log"
    benchmark:
        "out/{tool}{extra}_{b_id}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        ext_a="bam/bed/gff/vcf",
        tool="bedtools/coverage"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools coverage -a {input.a} -b {input.b} {params.extra} 2> {log}"
