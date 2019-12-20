rule bbmap_reformat:
    """
    Created:
        2019-03-21 11:03:53
    Aim:
        First used to select subpopulations from bam.
    Note:
        Currently does not work because input bam is treated as unpaired...
    Test:
        out/bbmap/reformat_nuc-length/samtools/sort/samtools/merge_three_samples/samtools/merge_five_runs/samtools/view_sam_to_bam/bowtie2/pe_mm10/ln/pe_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run163_run167_run184_run187_run205/s1-MNase_Spm_WT_band1_s2-MNase_Spm_WT_band3_s3-MNase_Spm_WT_band6.bam
    """
    input:
        "out/{filler}"
    output:
        "out/{tool}{extra}/{filler}"
    log:
        "out/{tool}{extra}/{filler}.log"
    benchmark:
        "out/{tool}{extra}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "bbmap/reformat"
    conda:
        "../envs/bbmap.yaml"
    shell:
        "reformat.sh in={input} out={output} {params.extra} &> {log}"

rule bbmap_se:
    """
    Created:
        2019-12-02 01:30:30
    Aim:
        First used to select subpopulations from bam.
        Trying this:
        https://www.biostars.org/p/225198/
    Note:
    Test:
        out/bbmap/se_fa-genome-GRCh38-r94-chr1/ln/alias/sst/all_samples/fastq/Jurkat_SRR1057274_H3K27ac.sam.gz
    """
    input:
        fq="out/{filler}.{ext}",
        fa= lambda wildcards: eval(config['ids'][wildcards.fa_genome_id])
    output:
        sam="out/{tool}_{ext}{extra}_{fa_genome_id}/{filler}.sam.gz"
    log:
        "out/{tool}_{ext}{extra}_{fa_genome_id}/{filler}.log"
    benchmark:
        "out/{tool}_{ext}{extra}_{fa_genome_id}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    threads:
        MAX_THREADS
    wildcard_constraints:
        tool = "bbmap/se",
        ext = "fastq.gz|fasta" #ADD MORE HERE IF NEEDED
    conda:
        "../envs/bbmap.yaml"
    shell:
        "bbmap.sh ref={input.fa} in={input.fq} out={output.sam} {params.extra} t={threads} maxindel=100k &> {log}"


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
        fa= lambda wildcards: eval(config['ids'][wildcards.fa_genome_id])
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

rule bbmap_readlength:
    """
    Created:
        2019-12-02 01:30:30
    Aim:
        First used to select subpopulations from bam.
        Trying this:
        https://www.biostars.org/p/225198/
    Note:
    Test:
        out/bbmap/readlength_csfasta.gz/ln/alias/experiments/solid_rnaseq_tall/ALL-T-2.10_12_F5.txt
    """
    input:
        "out/{filler}.{ext}",
    output:
        txt="out/{tool}_{ext}{extra}/{filler}.txt"
    log:
        "out/{tool}_{ext}{extra}/{filler}.log"
    benchmark:
        "out/{tool}_{ext}{extra}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "bbmap/readlength",
        ext = "fa|fq|sam|fastq.gz|csfasta.gz" #ADD HERE OTHER SUPPORTED EXTENSIONS I CURRENTLY HAVE NO INTEREST IN
    conda:
        "../envs/bbmap.yaml"
    shell:
        "readlength.sh in={input} out={output.txt} {params.extra} &> {log}"

