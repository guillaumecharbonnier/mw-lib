rule sra_tools_fastq_dump_se_extra:
    """
    Created:
        2018-03-19 16:54:48
    Aim:
        I have just realized that my url in src/snakemake/table/spermiogenesis_samples.tsv are broken.
        But fastq-dump is able to directly download sra given an id so let's try it!
    Test:
        out/sra-tools/fastq-dump_se/SRR1202037.fastq
        out/sra-tools/fastq-dump_se/SRR2096391.fastq
        out/sra-tools/fastq-dump_se/SRR646457.fastq.gz
    """
    output:
        fastq="out/{tool}{extra}/{srr}.{ext}"
    params:
        #outdir="out/{tool}{extra}/",
        extra = params_extra,
        gz = lambda wildcards: "--gzip" if wildcards.ext == 'fastq.gz' else ""
    wildcard_constraints:
        tool="sra-tools/fastq-dump_se",
        srr="SRR[0-9]+",
        ext="fastq|fastq.gz"
    conda:
        "../envs/sra-tools.yaml"
    threads:
        1
    shell:
        "fastq-dump {params} --outdir `dirname {output.fastq}` {wildcards.srr}"



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
    params:
        extra= params_extra,
        gz = lambda wildcards: "--gzip" if wildcards.ext == 'fastq.gz' else ""
    wildcard_constraints:
        tool="sra-tools/fastq-dump_pe",
        srr="SRR[0-9]+",
        ext="fastq|fastq.gz"
    conda:
        "../envs/sra-tools.yaml"
    threads:
        1
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

