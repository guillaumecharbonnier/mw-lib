rule pbgzip_tabix:
    """
    Created:
        2020-01-14 17:54:18
    Aim:
        Genralization of rule pbgzip_tabix_bed. 
        2019-10-15 18:14:50
    Modfied:
        2020-01-14 17:53:10 - Generalized to all files to use for vcf
    Aim:
        Test the combination of pbgzip and tabix to compress and index bed files.
        Originally done to test if Cyverse can recognize automatically such bed compressed files as bedgz.

    Test:
        out/pbgzip_tabix/bcftools/mpileup_fa-genome-hg19-main-chr/samtools/index/samtools/sort/samtools/view_sam_to_bam/bowtie2/se_fasta_--rfg_1,1_-k_1_-q_-f_hg19/edena/assembling_-d_20_-c_20_-minCoverage_5/edena/overlapping/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_H3K27ac_contigs.vcf.gz.tbi

        out/pbgzip_tabix/platypus/callVariants_--assemble=1_--minPosterior=0_fa-genome-hg19-main-chr_vcfgz-hg19-edena-contigs-with-inserts-gt-3bp-in-T11C-H3K27ac/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem2_se_bwa-index-hg19-main-chr-and-contigs-with-inserts-from-T11C-H3K27ac/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_H3K27ac.vcf.gz
        
    """
    input:
        "out/{filler}.{ext}",
    output:
        gz  = "out/pbgzip_tabix/{filler}.{ext}.gz",
        tbi = "out/pbgzip_tabix/{filler}.{ext}.gz.tbi"
    params:
        tmp = "out/pbgzip_tabix/{filler}.{ext}"
    log:
              "out/pbgzip_tabix/{filler}.{ext}.log"
    benchmark:
              "out/pbgzip_tabix/{filler}.{ext}.benchmark.tsv"
    conda:
        "../envs/pbgzip_tabix.yaml"
    wildcard_constraints:
        ext="bed|sam|vcf|gff"
    threads:
        1
    shell:
        """
        (ln -srf {input} {params.tmp}
        pbgzip {params.tmp}
        tabix -p {wildcards.ext} {output.gz}
        ) &> {log}
        """

rule pbgzip_tabix_bed_legacy:
    """
    Created:
        2019-10-15 18:14:50
    Modfied:
        2020-01-14 17:53:10 - Added pbgzip_tabix rule which should put this one to legacy
    Aim:
        Test the combination of pbgzip and tabix to compress and index bed files.
        Originally done to test if Cyverse can recognize automatically such bed compressed files as bedgz.
    Test:

    """
    input:
        bed="out/{filler}.bed",
    output:
        bedgz="out/pbgzip_tabix_bed/{filler}.bed.gz",
        tbi="out/pbgzip_tabix_bed/{filler}.bed.gz.tbi"
    params:
        bed="out/pbgzip_tabix_bed/{filler}.bed"
    log:
            "out/pbgzip_tabix_bed/{filler}.log"
    benchmark:
            "out/pbgzip_tabix_bed/{filler}.benchmark.tsv"
    conda:
        "../envs/pbgzip_tabix.yaml"
    threads:
        1
    shell:
        """
        (ln -srf {input.bed} {params.bed}
        pbgzip {params.bed}
        tabix -p bed {output.bedgz}
        ) &> {log}
        """

