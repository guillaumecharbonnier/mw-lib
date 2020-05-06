rule bowtie2_build:
    """
    Modified:
        2017-05-06 14:53:46 - Tool is now installed with conda.
    Doc:
        http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
    Test:
    out/bowtie2-build/wget/ftp/ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz.1.bt2
    out/bowtie2-build/fa-genome-GRCh38-chr22.1.bt2
        #index_parts  = expand("out/bowtie2-build/{{filler}}.{parts}.bt2",
    """
    input:
        #fasta="out/{filler}",
        fasta = lambda wildcards: eval(config['ids'][wildcards.fa_genome_id])
    output:
        index_parts  = expand("out/bowtie2-build/{{fa_genome_id}}.{parts}.bt2",
            parts=["1","2","3","4","rev.1","rev.2"])
    #log:
    #                      "out/{tool}_{ext}{extra}_{index}/{filler}.log"
    params:
        index="out/bowtie2-build/{fa_genome_id}"
        #index="out/bowtie2-build/{filler}"
    conda:
        "../envs/bowtie2.yaml"
    threads:
        MAX_THREADS
    shell: 
        """
        bowtie2-build {input.fasta} {params.index}
        """


