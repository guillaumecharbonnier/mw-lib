

rule sickle_se_extra:
    """
    Created:
        2019-01-31 21:19:42
    Aim:
        Trim reads based on sequencing quality threshold for single-end reads.
    Note:
        Recent Illumina base-caller usually provides Sanger quality type.
    Test:
        out/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR1202037.fastq
        out/sickle/se/sra-tools/fastq-dump_se/SRR1202037.fastq
    """
    input:
        fastq = "out/{filler}.{ext}"
    output:
        fastq = "out/{tool}{extra}/{filler}.{ext}"
    params:
        extra = params_extra,
        gz = lambda wildcards: "-g" if wildcards.ext == 'fastq.gz' else ""
    log:
        "out/{tool}{extra}/{filler}.{ext}.log"
    wildcard_constraints:
        tool="sickle/se",
        ext="fastq|fastq.gz"
    conda:
        "../envs/sickle.yaml"
    threads:
        1
    shell:
        "sickle se -f {input.fastq} {params} -o {output.fastq} &> {log}"

rule sickle_pe_extra:
    """
    Created:
        2019-01-31 21:34:32
    Aim:
        Trim reads based on sequencing quality threshold for paired-end reads.
    Test:
        out/sickle/pe/sra-tools/fastq-dump_pe/SRR3938832_2.fastq
        out/sickle/pe_-t_sanger_-q_30/sra-tools/fastq-dump_pe/SRR3938832_2.fastq
    """
    input:
        m1="out/{filler}_1.{ext}",
        m2="out/{filler}_2.{ext}"
    output:
        m1=    "out/{tool}{extra}/{filler}_1.{ext}",
        m2=    "out/{tool}{extra}/{filler}_2.{ext}",
        single="out/{tool}{extra}/{filler}_single.{ext}"
    params:
        extra = params_extra,
        gz = lambda wildcards: "-g" if wildcards.ext == 'fastq.gz' else ""
    log:
        "out/{tool}{extra}/{filler}.{ext}.log"
    wildcard_constraints:
        tool="sickle/pe",
        ext="fastq|fastq.gz"
    conda:
        "../envs/sickle.yaml"
    shell:
        "sickle pe -f {input.m1} -r {input.m2} "
        "-o {output.m1} -p {output.m2} "
        "-s {output.single} {params} &> {log}"

