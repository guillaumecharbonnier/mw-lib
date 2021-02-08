rule pblat:
    """
    Created:
        2019-11-28 15:29:59
    Doc:
        pblat -threads=<int> [other options] <database.fasta> <query.fasta> <out.psl>
    Note:
        pblat segfault on Sacapus node if database is the whole hg19.2bit file.
    Test:
        #out/pblat/dev_2bit-hg19-chr1/awk/sam_to_fasta/samtools/view_bam_to_sam/bedtools/intersect_-v_-b_bed-hg19-refgene-exons/out/samtools/index/samtools/sort/samtools/view_sam_to_bam/awk/extract_reads_with_insertions/bowtie2/se_-k_1_-q_hg19/bowtie/se_--chunkmbs_256_--best_--strata_-m_1_-n_2_ebwt-hg19/sickle/se_-t_sanger_-q_20/ln/alias/sst/all_samples/fastq/TH134_CD34_H3K27ac_unmapped.psl
        out/pblat/_-minScore=0_-stepSize=1_2bit-hg19-chr1/awk/sam_to_fasta/samtools/view_bam_to_sam/bedtools/intersect_-v_-b_bed-hg19-refgene-exons/out/samtools/index/samtools/sort/samtools/view_sam_to_bam/awk/extract_reads_with_insertions/bowtie2/se_-k_1_-q_hg19/bowtie/se_--chunkmbs_256_--best_--strata_-m_1_-n_2_ebwt-hg19/ln/alias/sst/all_samples/fastq/TH134_CD34_H3K27ac_unmapped.psl
        out/pblat/_-minScore=0_-stepSize=1_2bit-hg19-chr1/edena/assembling_-d_20_-c_20_-minCoverage_5/edena/overlapping/gunzip/to-stdout/bowtie/se_--chunkmbs_256_--best_--strata_-m_1_-n_2_ebwt-hg19-main-chr/ln/alias/sst/all_samples/fastq/T11C_H3K27ac_unmapped_contigs.psl
    """
    input:
        fasta="out/{filler}.fasta",
        ref=lambda wildcards: eval(mwconf['ids'][wildcards.ref_id])
    output:
        psl="out/{tool}{extra}_{ref_id}/{filler}.psl"
    log:
            "out/{tool}{extra}_{ref_id}/{filler}.log"
    benchmark:
            "out/{tool}{extra}_{ref_id}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="pblat/"
    conda:
        "../../../../mw-lib/src/snakemake/envs/pblat.yaml"
    threads:
        MAX_THREADS
    shell:
        "pblat -threads={threads} {params.extra} {input.ref} {input.fasta} {output.psl} &> {log}"
        #"pblat -threads={threads} -minScore=0 -stepSize=1 {input.ref} {input.fasta} {output.psl}"

