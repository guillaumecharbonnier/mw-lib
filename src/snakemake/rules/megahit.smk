rule megahit_se:
    """
    Created:
        2020-08-23 10:49:12
    Aim:
        Try an alternative to Edena when it fails to run because the input file is to big, e.g. T11C_all_DNA_samples without filtering.
    Doc:
        http://cab.spbu.ru/files/release3.14.0/manual.html
    Test:
        out/megahit/se/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac/final.contigs.fa
    """
    input:
        fq="out/{filler}.fastq.gz"
    output:
        expand("out/{{tool}}{{extra}}/{{filler}}/{filename}",filename=["final.contigs.fa"]) #Many more files here
    log:
        "out/{tool}{extra}/{filler}/log"
    benchmark:
        "out/{tool}{extra}/{filler}/benchmark.tsv"
    params:
        outdir="out/{tool}{extra}/{filler}",
        extra = params_extra
    wildcard_constraints:
        tool="megahit/se"
    conda:
        "../envs/megahit.yaml"
    threads:
        1 # Spades scales memory consumption with number of threads, and can use huge amount of RAM.
    shell:
        """
        rm -rf {params.outdir} #Megahit does not run if outdir exists
        megahit -t {threads} {params.extra} -r {input.fq} -o {params.outdir} #&> {log}
        """

rule megahit_pe:
    """
    Created:
        2020-08-17 09:44:11
    Doc:
        http://cab.spbu.ru/files/release3.14.0/manual.html
    Note:
        If memory consumption is too high, fastq may be filtered using bbnorm first.
        http://seqanswers.com/forums/showthread.php?t=49763
    Test:
        out/megahit/pe/ln/alias/sst/all_samples/fastq/856_H3K27ac/final.contigs.fa
    """
    input:
        fq_1="out/{filler}_1.fastq.gz",
        fq_2="out/{filler}_2.fastq.gz",
    output:
        outputs=expand("out/{{tool}}{{extra}}/{{filler}}/{filename}",filename=["final.contigs.fa"]) #Many more files here
        #expand("out/{{tool}}{{extra}}/{{filler}}/{filename}",filename=["contigs.fasta", "scaffolds.fasta"]) #Many more files here
        #expand("out/{{tool}}{{extra}}/{{filler}}/{filename}",filename=["transcripts.fasta"]) # Many more files here
    log:
        true="out/{tool}{extra}/{filler}/log"
    benchmark:
        "out/{tool}{extra}/{filler}/benchmark.tsv"
    params:
        outdir = "out/{tool}{extra}/{filler}",
        extra = params_extra
    wildcard_constraints:
        tool="megahit/pe"
    conda:
        "../envs/megahit.yaml"
    threads:
        1 # Spades scales memory consumption with number of threads, and can use huge amount of RAM.
    shell:
        """
        rm -rf {params.outdir} #Megahit does not run if outdir exists
        megahit -t {threads} {params.extra} \
            -1 {input.fq_1} -2 {input.fq_2} \
            -o {params.outdir} # &> {log.true}
        """

