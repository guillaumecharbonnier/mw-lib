rule bedtools_getfasta_extra:
    """
    Modified:
        2017-02-06 15:05:57 - Modified input and output to allow the use of mm9 assembly and other beds. Need to write a rule to download mm9 fasta.
    Doc:
        https://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html
    Test:
        out/bedtools/getfasta_fa-genome-mm10/macs2/callpeak_--broad/samtools/index/samtools/sort/samtools/view_sam_to_bam_-q_30/bowtie2/se_mm10/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR3126243_over_SRR3126242_peaks.fa
    """
    input:
        #fasta = input_fa_genome,
        fasta = lambda wildcards: eval(config['ids'][wildcards.fa_genome_id]),
        #fai   = lambda wildcards: [path + '.fai' for path in eval(config['ids'][wildcards.bam_list_id])],
        bed = "out/{filler}.bed"
    output:
        fasta="out/{tool}{extra}_{fa_genome_id}/{filler}.fa"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="bedtools/getfasta"
    conda:
        "../envs/bedtools.yaml"
    threads:
        1
    shell:
        "bedtools getfasta -fi {input.fasta} -bed {input.bed} -fo {output.fasta} {params.extra}"
