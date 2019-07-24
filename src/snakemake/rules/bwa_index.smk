rule bwa_index:
    """
    Created:
        2017-05-14 13:14:45
    Doc:
        http://bio-bwa.sourceforge.net/bwa.shtml
    Test:
        out/bwa/index/fa-genome-GRCh38-Blueprint.amb
    """
    input:
        fasta = lambda wildcards: eval(config['ids'][wildcards.fa_genome_id])
        #fai   = lambda wildcards: [path + '.fai' for path in eval(config['ids'][wildcards.bam_list_id])],
        #genome=input_fa_genome
    output:
        index=expand("out/bwa/index/{{fa_genome_id}}.{ext}", ext=["amb","ann","bwt","pac","sa"])
    params:
        prefix="out/bwa/index/{fa_genome_id}"
    conda:
        "../envs/bwa.yaml"
    threads:
        1
    shell:
        "bwa index -p {params.prefix} {input.fasta}"

