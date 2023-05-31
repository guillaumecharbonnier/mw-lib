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
    Note:
        This rule was perfectly fine unless run on Bigmemorix on the sacapus mount point:
        leading to errors:
        Output directory 'out/fastq_screen/filter/sickle/se_-t_sanger_-q_20/ln/alias/sst/all_samples/fastq' is not writable, please adjust configuration.
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
        extra = params_extra,
        # fastq_screen_config="fastq_screen.conf",
        fastq_screen_config="out/sed/edit_fastq_screen_conf/fastq_screen/get_genomes/FastQ_Screen_Genomes/fastq_screen.conf",
        subset=100000,
        aligner='bowtie2'
    wildcard_constraints:
        tool = "fastq_screen/filter"
    threads:
        MAX_THREADS
    shell:
        # rm needed because fastq_screen send an error message 
        # if the output folder already exists from a previous run
        #"rm -rf `dirname {output.txt}`; "
        "fastq_screen {params.extra} "
        "--conf {input.conf} --threads {threads} "
        "--outdir `dirname {output.txt}` {input.fastq} &> {log}"

rule fastq_screen_tests:
    input:
        "out/fastq_screen/filter/sickle/se_-t_sanger_-q_20/ln/alias/sst/all_samples/fastq/TH134_CD34_H3K27ac_screen.txt"

rule fastq_screen_wrapper:
    input:
        fastq = "out/{filler}.fastq.gz",
        conf = "out/sed/edit_fastq_screen_conf/fastq_screen/get_genomes/FastQ_Screen_Genomes/fastq_screen.conf"
    output:
        txt="out/{tool}{extra}/{filler}_screen.txt",
        png="out/{tool}{extra}/{filler}_screen.png",
        # html="out/{tool}{extra}/{filler}_screen.html"
        #html is not exported in current wrapper version
    log:
        "out/{tool}{extra}/{filler}/log"
    benchmark:
        "out/{tool}{extra}/{filler}/benchmark.tsv"
    params:
        extra = params_extra,
        # fastq_screen_config="fastq_screen.conf",
        fastq_screen_config="out/sed/edit_fastq_screen_conf/fastq_screen/get_genomes/FastQ_Screen_Genomes/fastq_screen.conf",
        subset=100000,
        aligner='bowtie2'
    wildcard_constraints:
        tool = "fastq_screen/wrapper"
    # input:
    #     "samples/{sample}.fastq"
    # output:
    #     txt="qc/{sample}.fastq_screen.txt",
    #     png="qc/{sample}.fastq_screen.png"
    threads: 8
    wrapper:
        "v1.28.0/bio/fastq_screen"
