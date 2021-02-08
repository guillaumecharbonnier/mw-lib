rule bbmap_callvariants:
    """
    Created:
        2019-12-02 01:30:30
    Aim:
        First used to select subpopulations from bam.
        Trying this:
        https://www.biostars.org/p/225198/
    Note:
    Test:
        out/bbmap/callvariants_fa-genome-GRCh38-r94-chr1/bbmap/se_fa-genome-GRCh38-r94-chr1/ln/alias/sst/all_samples/fastq/Jurkat_SRR1057274_H3K27ac.vcf
        out/bbmap/callvariants_fa-genome-hg19-main-chr/bbmap/se_fasta_fa-genome-hg19-main-chr/edena/assembling_-d_20_-c_20_-minCoverage_5/edena/overlapping/gunzip/to-stdout/bowtie/se_--chunkmbs_256_--best_--strata_-m_1_-n_2_ebwt-hg19-main-chr/ln/alias/sst/all_samples/fastq/T11C_H3K27ac_unmapped_contigs.vcf
    """
    input:
        sam="out/{filler}.sam.gz",
        fa= lambda wildcards: eval(mwconf['ids'][wildcards.fa_genome_id])
    output:
        vcf="out/{tool}{extra}_{fa_genome_id}/{filler}.vcf"
    log:
        "out/{tool}{extra}_{fa_genome_id}/{filler}.log"
    benchmark:
        "out/{tool}{extra}_{fa_genome_id}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    threads:
        MAX_THREADS
    wildcard_constraints:
        tool = "bbmap/callvariants"
    conda:
        "../envs/bbmap.yaml"
    shell:
        "callvariants.sh ref={input.fa} in={input.sam} vcf={output.vcf} {params.extra} t={threads} ploidy=2 &> {log}"

