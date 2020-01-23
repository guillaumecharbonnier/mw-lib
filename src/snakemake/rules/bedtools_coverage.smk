rule bedtools_coverage:
    """
    Created:
        2016-08-24 14h59
    Aim:
        The bedtools coverage tool computes both the depth and breadth of coverage of features in file B on the features in file A. 
        For example, bedtools coverage can compute the coverage of sequence alignments (file B) across 1 kilobase (arbitrary) windows (file A) tiling a genome of interest.
        One advantage that bedtools coverage offers is that it not only counts the number of features that overlap an interval in file A, it also computes the fraction of bases in the interval in A that were overlapped by one or more features.
        Thus, bedtools coverage also computes the breadth of coverage observed for each interval in A.
    Doc:
        https://bedtools.readthedocs.io/en/latest/content/tools/coverage.html
    Note:
        Bedtools doc mention usage like this:
        -b <FILE1, FILE2, ..., FILEN>
        but actually FILES should be separated by a space only.
    Test:
        out/bedtools/coverage_bam-hg19-T11C-H3K27ac-dev1/bcftools/mpileup_fa-genome-hg19-main-chr/samtools/index/samtools/sort/samtools/view_sam_to_bam/bowtie2/se_fasta_--rfg_1,1_-k_1_-q_-f_hg19/edena/assembling_-d_20_-c_20_-minCoverage_5/edena/overlapping/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_H3K27ac_contigs.vcf

        Debug: out/bcftools/mpileup_fa-genome-hg19-main-chr/samtools/index/samtools/sort/samtools/view_sam_to_bam/bowtie2/se_fasta_--rfg_1,1_-k_1_-q_-f_hg19/edena/assembling_-d_20_-c_20_-minCoverage_5/edena/overlapping/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_H3K27ac_contigs.vcf
    """
    input:
        a="out/{filler}.{ext_a}",
        b=lambda wildcards: eval(config['ids'][wildcards.b_id])
    output:
        "out/{tool}{extra}_{b_id}/{filler}.{ext_a}"
    log:
        "out/{tool}{extra}_{b_id}/{filler}.{ext_a}.log"
    benchmark:
        "out/{tool}{extra}_{b_id}/{filler}.{ext_a}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        ext_a="bam|bed|gff|vcf",
        tool="bedtools/coverage"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools coverage -a {input.a} -b {input.b} {params.extra} 1> {output} 2> {log}"


rule bedtools_coverage_test:
    """
    Created:
        2016-08-24 14h59
    Aim:
        The bedtools coverage tool computes both the depth and breadth of coverage of features in file B on the features in file A. 
        For example, bedtools coverage can compute the coverage of sequence alignments (file B) across 1 kilobase (arbitrary) windows (file A) tiling a genome of interest.
        One advantage that bedtools coverage offers is that it not only counts the number of features that overlap an interval in file A, it also computes the fraction of bases in the interval in A that were overlapped by one or more features.
        Thus, bedtools coverage also computes the breadth of coverage observed for each interval in A.
    Doc:
        https://bedtools.readthedocs.io/en/latest/content/tools/coverage.html
    Note:
        Draft, maybe deprecated for my current need in comparison with deepTools multbigwigsummary + plotCorrelation
    Test:
        out/bedtools/coverage_test_bam-hg19-T11C-H3K27ac-dev1/bcftools/mpileup_fa-genome-hg19-main-chr/samtools/index/samtools/sort/samtools/view_sam_to_bam/bowtie2/se_fasta_--rfg_1,1_-k_1_-q_-f_hg19/edena/assembling_-d_20_-c_20_-minCoverage_5/edena/overlapping/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_H3K27ac_contigs.vcf

        Debug: out/bcftools/mpileup_fa-genome-hg19-main-chr/samtools/index/samtools/sort/samtools/view_sam_to_bam/bowtie2/se_fasta_--rfg_1,1_-k_1_-q_-f_hg19/edena/assembling_-d_20_-c_20_-minCoverage_5/edena/overlapping/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_H3K27ac_contigs.vcf
    """
    input:
        a="out/{filler}.{ext_a}",
        b=lambda wildcards: eval(config['ids'][wildcards.b_id])
    output:
        "out/{tool}{extra}_{b_id}/{filler}.{ext_a}"
    log:
        "out/{tool}{extra}_{b_id}/{filler}.{ext_a}.log"
    benchmark:
        "out/{tool}{extra}_{b_id}/{filler}.{ext_a}.benchmark.tsv"
    params:
        extra = params_extra,
    wildcard_constraints:
        ext_a="bam|bed|gff|vcf",
        tool="bedtools/test_coverage"
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        INPUT_B=`echo {input.b} | sed 's/ /, /g'`
        bedtools coverage -a {input.a} -b $INPUT_B {params.extra} 1> {output} 2> {log}
        """
