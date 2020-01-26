rule freebayes_fa_bam:
    """
    Created:
        2020-01-22 22:18:19
    Doc:
        https://github.com/ekg/freebayes
    Note:
    Test:
        out/pbgzip_tabix/freebayes/_fa-genome-hg19-main-chr/samtools/index/abra2/_--single_fa-genome-hg19-main-chr/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_se_fa-genome-hg19-main-chr/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_H3K27ac.vcf.gz
    """
    input:
        bam="out/{filler}.bam",
        fa= lambda wildcards: eval(config['ids'][wildcards.fa_genome_id])
    output:
        vcf="out/{tool}{extra}_{fa_genome_id}/{filler}.vcf"
    log:
            "out/{tool}{extra}_{fa_genome_id}/{filler}.log"
    benchmark:
            "out/{tool}{extra}_{fa_genome_id}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    conda:
        "../envs/freebayes.yaml"
    wildcard_constraints:
        tool="freebayes/"
    threads:
        1
    shell:
        """
        freebayes {params.extra} -f {input.fa} {input.bam} 1> {output.vcf} 2> {log}
        """

rule freebayes_fa_vcf_bam:
    """
    Created:
        2020-01-22 14:54:54
    Doc:
        https://github.com/mozack/abra2
    Note:
        Single end data must have --single argument. Paired end is assumed by default.
    Test:
        out/pbgzip_tabix/freebayes/_fa-genome-hg19-main-chr_vcfgz-hg19-edena-contigs-with-indels-gt-3bp-in-T11C-H3K27ac/samtools/index/abra2/_--single_fa-genome-hg19-main-chr_vcf-hg19-edena-contigs-with-indels-gt-3bp-in-T11C-H3K27ac/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_se_fa-genome-hg19-main-chr/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_H3K27ac.vcf.gz
    """
    input:
        bam="out/{filler}.bam",
        fa  = lambda wildcards: eval(config['ids'][wildcards.fa_genome_id]),
        vcfgz = lambda wildcards: eval(config['ids'][wildcards.vcfgz_id])
    output:
        vcf="out/{tool}{extra}_{fa_genome_id}_{vcfgz_id}/{filler}.vcf"
    log:
            "out/{tool}{extra}_{fa_genome_id}_{vcfgz_id}/{filler}.log"
    benchmark:
            "out/{tool}{extra}_{fa_genome_id}_{vcfgz_id}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    conda:
        "../envs/freebayes.yaml"
    wildcard_constraints:
        tool="freebayes/"
    threads:
        1
    shell:
        """
        freebayes {params.extra} -f {input.fa} -@ {input.vcfgz} {input.bam} 1> {output.vcf} 2> {log}
        """


