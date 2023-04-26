rule fastq_screen_get_genomes:
    """
    Created:
        2019-03-06 13:18:42
    Aim:
        Obtaining reference genomes
    Doc:
        https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/_build/html/index.html#obtaining-reference-genomes
    """
    output:
        "out/fastq_screen/get_genomes/FastQ_Screen_Genomes/fastq_screen.conf"
    conda:
        "../envs/fastq_screen.yaml"
    shell:
        "rm -rf out/fastq_screen/get_genomes/FastQ_Screen_Genomes; fastq_screen --get_genomes --outdir out/fastq_screen/get_genomes"

rule sed_edit_fastq_screen_conf:
    """
    """
    input:
        "out/fastq_screen/get_genomes/FastQ_Screen_Genomes/fastq_screen.conf"
    output:
        "out/sed/edit_fastq_screen_conf/fastq_screen/get_genomes/FastQ_Screen_Genomes/fastq_screen.conf"
    params:
        genome_path = "out/fastq_screen/get_genomes"
    conda:
        "../envs/coreutils.yaml"
    shell:
        "sed -e 's/THREADS		7/THREADS		16/' "
        "-e 's|/bi/scratch/wingetts/FASTQ_Screen_Paper/FastQ_Screen_Genomes|{params.genome_path}|' {input} > {output}"

rule fastq_screen_filter:
    """
    Created:
        2019-03-06 13:18:42
    Aim:
        Obtaining reference genomes
    Doc:
        https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/_build/html/index.html
    Test:
        out/fastq_screen/filter/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR3126243_screen.txt
    """
    input:
        fastq = "out/{filler}.fastq.gz",
        conf = "out/sed/edit_fastq_screen_conf/fastq_screen/get_genomes/FastQ_Screen_Genomes/fastq_screen.conf"
    output:
        txt="out/{tool}{extra}/{filler}_screen.txt",
        png="out/{tool}{extra}/{filler}_screen.png",
        html="out/{tool}{extra}/{filler}_screen.html"
    log:
        "out/{tool}{extra}/{filler}/log"
    benchmark:
        "out/{tool}{extra}/{filler}/benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "fastq_screen/filter"
    conda:
        "../envs/fastq_screen.yaml"
    threads:
        MAX_THREADS
    shell:
        # rm needed because fastq_screen send an error message 
        # if the output folder already exists from a previous run
        #"rm -rf `dirname {output.txt}`; "
        "fastq_screen {params.extra} "
        "--conf {input.conf} --threads {threads} "
        "--outdir `dirname {output.txt}` {input.fastq} &> {log}"


