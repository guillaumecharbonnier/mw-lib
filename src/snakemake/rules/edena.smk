rule edena_overlapping:
    """
    Modified: 
        2020-08-17 10:05:36 Added '_se' in tool name. Should lead to deprecation in old code. TODO: Add '_se' everywhere.
    Doc:
        https://oit.ua.edu/wp-content/uploads/2016/10/edena_referencemanual120926.pdf
    Note:
        Reads provided to edena should have same length, hence quality trimming, e.g. with sickle, should be skipped
    Test:
        out/edena/overlapping_se/awk/sam_to_fastq/samtools/view_bam_to_sam/bedtools/intersect_-v_-b_bed-hg19-refgene-exons/samtools/index/samtools/sort/samtools/view_sam_to_bam/awk/extract_reads_with_insertions/bowtie2/se_-k_1_-q_hg19/bowtie/se_--chunkmbs_256_--best_--strata_-m_1_-n_2_ebwt-hg19/ln/alias/sst/all_samples/fastq/TH134_CD34_H3K27ac_unmapped.ovl
    """
    input:
        fastq="out/{filler}.fastq"
    output:
        "out/{tool}{extra}/{filler}.ovl"
    log:
        "out/{tool}{extra}/{filler}.log"
    benchmark:
        "out/{tool}{extra}/{filler}.benchmark.tsv"
    params:
        outdir="out/{tool}{extra}/{filler}",
        extra = params_extra
    wildcard_constraints:
        tool="edena/overlapping_se"
    threads:
        1
        # Actually only use one thread even if more are provided
        #MAX_THREADS
    conda:
        "../envs/edena.yaml"
    shell:
        "edena -nThreads {threads} {params.extra} -r {input.fastq} -p {params.outdir} &> {log}"


rule edena_overlapping_pe:
    """
    Created:
        2020-08-17 09:59:35
    Doc:
        https://oit.ua.edu/wp-content/uploads/2016/10/edena_referencemanual120926.pdf
    Note:
        Reads provided to edena should have same length, hence quality trimming, e.g. with sickle, should be skipped
    Test:
        out/edena/overlapping_pe/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/856_H3K27ac.ovl
    """
    input:
        fq_1="out/{filler}_1.fastq",
        fq_2="out/{filler}_2.fastq"
    output:
        "out/{tool}{extra}/{filler}.ovl"
    log:
        "out/{tool}{extra}/{filler}.log"
    benchmark:
        "out/{tool}{extra}/{filler}.benchmark.tsv"
    params:
        outdir="out/{tool}{extra}/{filler}",
        extra = params_extra
    wildcard_constraints:
        tool="edena/overlapping_pe"
    threads:
        1
        # Actually only use one thread even if more are provided
        #MAX_THREADS
    conda:
        "../envs/edena.yaml"
    shell:
        "edena -nThreads {threads} {params.extra} -DRpairs {input.fq_1} {input.fq_2} -p {params.outdir} &> {log}"

rule edena_assembling:
    """
    Doc:
        https://oit.ua.edu/wp-content/uploads/2016/10/edena_referencemanual120926.pdf
    Note:
        Reads provided to edena should have same length, hence quality trimming, e.g. with sickle, should be skipped
    Test:
        out/edena/assembling_-d_20_-c_20_-minCoverage_5/edena/overlapping/awk/sam_to_fastq/samtools/view_bam_to_sam/bedtools/intersect_-v_-b_bed-hg19-refgene-exons/samtools/index/samtools/sort/samtools/view_sam_to_bam/awk/extract_reads_with_insertions/bowtie2/se_-k_1_-q_hg19/bowtie/se_--chunkmbs_256_--best_--strata_-m_1_-n_2_ebwt-hg19/ln/alias/sst/all_samples/fastq/TH134_CD34_H3K27ac_unmapped_contigs.fasta

        out/edena/assembling/edena/overlapping/gunzip/to-stdout/bowtie/se_--chunkmbs_256_--best_--strata_-m_1_-n_2_ebwt-hg19/ln/alias/sst/all_samples/fastq/Jurkat_SRR1057274_H3K27ac_unmapped_contigs.fasta
        out/edena/assembling/edena/overlapping/gunzip/to-stdout/bowtie/se_--chunkmbs_256_--best_--strata_-m_1_-n_2_ebwt-hg19-main-chr/ln/alias/sst/all_samples/fastq/Jurkat_SRR1057274_H3K27ac_unmapped_contigs.fasta
        out/edena/assembling_-d_20_-c_20_-minCoverage_5/edena/overlapping/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/Jurkat_SRR1057274_H3K27ac_contigs.fasta

    """
    input:
        ovl="out/{filler}.ovl"
    output:
        "out/{tool}{extra}/{filler}_contigs.fasta",
        "out/{tool}{extra}/{filler}_contigs.lay",
    log:
        "out/{tool}{extra}/{filler}.log"
    benchmark:
        "out/{tool}{extra}/{filler}.benchmark.tsv"
    params:
        outdir="out/{tool}{extra}/{filler}",
        extra = params_extra
    wildcard_constraints:
        tool="edena/assembling"
    threads:
        MAX_THREADS
    conda:
        "../envs/edena.yaml"
    shell:
        "edena -nThreads {threads} {params.extra} -e {input.ovl} -p {params.outdir}"

# After assembling:
# /opt/bao/bin/bowtie2/bowtie2  --rfg 1,1 -p 5 -k 1 -q -f -x /opt/bao/bin/bowtie2/indexes/hg19/hg19 -U $file.unmapped/$file.unmapped_contigs.fasta -S $file.unmapped/$file.unmapped_contigs.sam
# out/bowtie2/se_--rfg_1,1_-k_1_-q_-f_hg19/edena/assembling_-d_20_-c_20_-minCoverage_5/edena/overlapping/awk/sam_to_fastq/samtools/view_bam_to_sam/bedtools/intersect_-v_-b_bed-hg19-refgene-exons/samtools/index/samtools/sort/samtools/view_sam_to_bam/awk/extract_reads_with_insertions/bowtie2/se_-k_1_-q_hg19/bowtie/se_--chunkmbs_256_--best_--strata_-m_1_-n_2_ebwt-hg19/ln/alias/sst/all_samples/fastq/TH134_CD34_H3K27ac_unmapped_contigs.sam




