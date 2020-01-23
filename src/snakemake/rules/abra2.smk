rule abra2_fa_bam:
    """
    Created:
        2020-01-22 14:54:54
    Doc:
        https://github.com/mozack/abra2
    Note:
        Single end data must have --single argument. Paired end is assumed by default.
    Test:
        out/abra2/_--single_fa-genome-hg19-main-chr/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_se_fa-genome-hg19-main-chr/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_H3K27ac.bam
    """
    input:
        bam="out/{filler}.bam",
        fa= lambda wildcards: eval(config['ids'][wildcards.fa_genome_id])
    output:
        bam="out/{tool}{extra}_{fa_genome_id}/{filler}.bam"
    log:
            "out/{tool}{extra}_{fa_genome_id}/{filler}.log"
    benchmark:
            "out/{tool}{extra}_{fa_genome_id}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    conda:
        "../envs/abra2.yaml"
    wildcard_constraints:
        tool="abra2/"
    threads:
        MAX_THREADS
    shell:
        """
        TMPDIR="`dirname {output}`/TMPDIR"
        echo $TMPDIR
        mkdir -p $TMPDIR

        # Check this regularly to know if export can be removed:
        # https://github.com/mozack/abra2/issues/25
        export LC_ALL=en_US.UTF-8
        abra2 {params.extra} --in {input.bam} --out {output.bam} --ref {input.fa} --threads {threads} --tmpdir $TMPDIR &> {log}
        rm -rf $TMPDIR
        """

rule abra2_fa_vcf_bam:
    """
    Created:
        2020-01-22 14:54:54
    Doc:
        https://github.com/mozack/abra2
    Note:
        Single end data must have --single argument. Paired end is assumed by default.
    Test:
        out/samtools/index/abra2/_--single_fa-genome-hg19-main-chr_vcf-hg19-edena-contigs-with-indels-gt-3bp-in-T11C-H3K27ac/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_se_fa-genome-hg19-main-chr/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_H3K27ac.bam
    """
    input:
        bam="out/{filler}.bam",
        fa  = lambda wildcards: eval(config['ids'][wildcards.fa_genome_id]),
        vcf = lambda wildcards: eval(config['ids'][wildcards.vcf_id])
    output:
        bam="out/{tool}{extra}_{fa_genome_id}_{vcf_id}/{filler}.bam"
    log:
            "out/{tool}{extra}_{fa_genome_id}_{vcf_id}/{filler}.log"
    benchmark:
            "out/{tool}{extra}_{fa_genome_id}_{vcf_id}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    conda:
        "../envs/abra2.yaml"
    wildcard_constraints:
        tool="abra2/"
    threads:
        MAX_THREADS
    shell:
        """
        TMPDIR="`dirname {output}`/TMPDIR"
        echo $TMPDIR
        mkdir -p $TMPDIR

        # Check this regularly to know if export can be removed:
        # https://github.com/mozack/abra2/issues/25
        export LC_ALL=en_US.UTF-8
        abra2 {params.extra} --in {input.bam} --in-vcf {input.vcf} --out {output.bam} --ref {input.fa} --threads {threads} --tmpdir $TMPDIR &> {log}
        rm -rf $TMPDIR
        """


