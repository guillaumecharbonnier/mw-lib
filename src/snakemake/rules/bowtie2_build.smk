rule bowtie2_build:
    """
    Modified:
        2017-05-06 14:53:46 - Tool is now installed with conda.
    Note:
        This rule will fail if fasta contains more than 4-billion nucleotides.
        See rule below instead.
    Doc:
        http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
    Test:
        out/bowtie2-build/wget/ftp/ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz.1.bt2
        out/bowtie2-build/fa-genome-GRCh38-chr22.1.bt2
        out/bowtie2-build/fa-ncrna-GRCh38-ensembl-r101.1.bt2
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
        bowtie2-build --threads {threads} {input.fasta} {params.index}
        """

rule bowtie2_build_large_index:
    """
    Created:
        2020-08-18 10:07:47
    Note:
        This variant of above rule is for fasta with more than 4-billion nucleotides when bowtie2-build output bt2l index instead of bt2 one.
    Doc:
        http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
    Test:
    """
    input:
        #fasta="out/{filler}",
        fasta = lambda wildcards: eval(config['ids'][wildcards.fa_genome_id])
    output:
        index_parts  = expand("out/bowtie2-build/{{fa_genome_id}}.{parts}.bt2l",
            parts=["1","2","3","4","rev.1","rev.2"])
    params:
        index="out/bowtie2-build/{fa_genome_id}"
    conda:
        "../envs/bowtie2.yaml"
    threads:
        MAX_THREADS
    shell: 
        """
        bowtie2-build --threads {threads} --large-index {input.fasta} {params.index}
        """


