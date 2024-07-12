rule sra_tools_fasterq_dump_se:
    """
    Test:
        out/sra-tools/fasterq-dump_se/GSM4566060.fastq.gz
        out/sra-tools/fasterq-dump_se/SRR1202037.fastq.gz
        out/sra-tools/fasterq-dump_se/ERR2308481.fastq.gz
    """
    output:
        fastq="out/sra-tools/fasterq-dump_se/{id}.fastq.gz"
    log:
        "out/sra-tools/fasterq-dump_se/{id}.log"
    benchmark:
        "out/sra-tools/fasterq-dump_se/{id}.benchmark.tsv"
    wildcard_constraints:
        id="(GSM|SRR|ERR)[0-9]+"
    conda:
        # "../envs/pysradb_sra-tools.yaml"
        "../envs/sra-tools.yaml"
    threads:
        1
    shell:
        """
        (
        DOWNLOAD_DIR="out/sra-tools/fasterq-dump_se/{wildcards.id}"
        # Maybe we can avoid rm here to allow resuming prefetch step if network issue occurs?
        #rm -rf $DOWNLOAD_DIR
        mkdir -p $DOWNLOAD_DIR
        cd $DOWNLOAD_DIR
        prefetch {wildcards.id}
        fasterq-dump {wildcards.id}
        gzip -c *.fastq > ../{wildcards.id}.fastq.gz
        rm *.fastq
        # cd ../
        # rm -rf {wildcards.id}
        ) &> {log}
        """

rule sra_tools_fasterq_dump_pe:
    """
    Test:
        out/sra-tools/fasterq-dump_pe/GSM4421303_2.fastq.gz
    """
    output:
        R1 = "out/sra-tools/fasterq-dump_pe/{id}_1.fastq.gz",
        R2 = "out/sra-tools/fasterq-dump_pe/{id}_2.fastq.gz"
    log:
        "out/sra-tools/fasterq-dump_pe/{id}.log"
    benchmark:
        "out/sra-tools/fasterq-dump_pe/{id}.benchmark.tsv"
    wildcard_constraints:
        id="(GSM|SRR|ERR)[0-9]+"
    conda:
        # "../envs/pysradb_sra-tools.yaml"
        "../envs/sra-tools.yaml"
    threads:
        1
    shell:
        """
        (
        DOWNLOAD_DIR="out/sra-tools/fasterq-dump_pe/{wildcards.id}"
        # Maybe we can avoid rm here to allow resuming prefetch step if network issue occurs?
        #rm -rf $DOWNLOAD_DIR
        mkdir -p $DOWNLOAD_DIR
        cd $DOWNLOAD_DIR
        prefetch {wildcards.id}
        fasterq-dump {wildcards.id}
        gzip -c *_1.fastq > ../{wildcards.id}_1.fastq.gz
        gzip -c *_2.fastq > ../{wildcards.id}_2.fastq.gz
        rm *.fastq
        # cd ../
        # rm -rf {wildcards.id}
        ) &> {log}
        """

rule sra_tools_fastq_dump_se_extra:
    """
    Created:
        2018-03-19 16:54:48
    Aim:
        I have just realized that my url in src/snakemake/table/spermiogenesis_samples.tsv are broken.
        But fastq-dump is able to directly download sra given an id so let's try it!
    Note:
        I have tested to query a SRX instead of a SRR but it does not work.
    Test:
        out/sra-tools/fastq-dump_se/SRR1202037.fastq
        out/sra-tools/fastq-dump_se/SRR2096391.fastq
        out/sra-tools/fastq-dump_se/SRR646457.fastq.gz
    """
    output:
        fastq="out/{tool}{extra}/{srr}.{ext}"
    log:
        "out/{tool}{extra}/{srr}.{ext}.log"
    benchmark:
        "out/{tool}{extra}/{srr}.{ext}.benchmark.tsv"
    params:
        #outdir="out/{tool}{extra}/",
        extra = params_extra,
        gz = lambda wildcards: "--gzip" if wildcards.ext == 'fastq.gz' else ""
    wildcard_constraints:
        tool="sra-tools/fastq-dump_se",
        srr="[ES]RR[0-9]+",
        ext="fastq|fastq.gz"
    conda:
        "../envs/sra-tools.yaml"
    resources:
        fastqdump_token=1
    threads:
        1
    shell:
        "fastq-dump {params} --outdir `dirname {output.fastq}` {wildcards.srr} &> {log}"

rule sra_tools_fastq_dump_pe_extra:
    """
    Created:
        2016-04-20 16h31
    Modified:
        2017-02-14 11:12:44 - Not tested
        2017-05-11 10:31:49 - Changed to legacy
        2017-05-11 10:47:18 - Use miniconda sra-tools now.
        2018-06-21 18:16:12 - Rely on srr and not an input sra.
    Note:
        --split-files allows to create one fastq file for each mate for paired-ends.
    Test:
        out/sra-tools/fastq-dump_pe/SRR3938832_2.fastq
        out/sra-tools/fastq-dump_pe/SRR1206271_2.fastq.gz
    """
    output:
        mates=expand("out/{{tool}}{{extra}}/{{srr}}_{mate}.{{ext}}", mate=["1","2"])
    log:
        "out/{tool}{extra}/{srr}.{ext}.log"
    benchmark:
        "out/{tool}{extra}/{srr}.{ext}.benchmark.tsv"
    params:
        extra= params_extra,
        gz = lambda wildcards: "--gzip" if wildcards.ext == 'fastq.gz' else ""
    wildcard_constraints:
        tool="sra-tools/fastq-dump_pe",
        srr="[SE]RR[0-9]+",
        #srr="SRR[0-9]+",
        ext="fastq|fastq.gz"
    conda:
        "../envs/sra-tools.yaml"
    threads:
        1
    resources:
        fastqdump_token=1
    shell:
        "fastq-dump --split-files {params} --outdir `dirname {output[0]}` {wildcards.srr}"

rule sra_tools_fastqdump_check_se_or_pe:
    """
    Created:
        2016-09-30 9h40
    Modified:
        2017-02-14 11:20:18
        2017-05-11 10:50:06 - Use miniconda sra-tools now.
    Aim:
        Print 1 if single-end and 2 if paired-end.
        Sometimes better than looking at <LIBRARY_LAYOUT> in this kind of XML when experiment merge both designs.
        https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=SRR1999050
    Test:
        out/sra-tools/fastq-dump_check_se_or_pe/SRR3938832.txt
        out/sra-tools/fastq-dump_check_se_or_pe/SRR1202037.txt
    """
    output:
        txt="out/sra-tools/fastq-dump_check_se_or_pe/{srr}.txt"
    conda:
        "../envs/sra-tools.yaml"
    shell:
        """
        fastq-dump -I -X 1 -Z --split-spot {wildcards.srr} 2>/dev/null \
            | awk '{{if(NR % 2 == 1) print substr($1,length($1),1)}}' \
            | uniq \
            | wc -l > {output.txt}
        """

