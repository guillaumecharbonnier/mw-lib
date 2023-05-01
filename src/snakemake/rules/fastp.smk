rule fastp_pe_extra:
    """
    Created:
        2023-04-26 11:08:56
    Aim:
        Trim reads based on sequencing quality threshold for paired-end reads.
    Test:
        out/fastp/pe/tidy_eurofins_filenames/K3_T1_NT_lib677680_10213_2_2.fastq.gz
        out/fastp/pe/sra-tools/fastq-dump_pe/SRR3938832_2.fastq
    """
    input:
        m1="out/{filler}_1.{ext}",
        m2="out/{filler}_2.{ext}"
    output:
        # single="out/{tool}{extra}/{filler}_single.{ext}"
        m1=    "out/{tool}{extra}/{filler}_1.{ext}",
        m2=    "out/{tool}{extra}/{filler}_2.{ext}"
    params:
        extra = params_extra,
    log:
        "out/{tool}{extra}/{filler}.{ext}.log"
    wildcard_constraints:
        tool="fastp/pe",
        ext="fastq|fastq.gz"
    conda:
        "../envs/fastp.yaml"
    shell:
        "fastp pe -i {input.m1} -I {input.m2} "
        "-o {output.m1} -O {output.m2} "
        "{params} &> {log}"



# rule fastp_se_extra:
#     """
#     Created:
#         2019-01-31 21:19:42
#     Aim:
#         Trim reads based on sequencing quality threshold for single-end reads.
#     Note:
#         Recent Illumina base-caller usually provides Sanger quality type.
#     Test:
#         out/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR1202037.fastq
#         out/sickle/se/sra-tools/fastq-dump_se/SRR1202037.fastq
#     """
#     input:
#         fastq = "out/{filler}.{ext}"
#     output:
#         fastq = "out/{tool}{extra}/{filler}.{ext}"
#     params:
#         extra = params_extra,
#         gz = lambda wildcards: "-g" if wildcards.ext == 'fastq.gz' else ""
#     log:
#         "out/{tool}{extra}/{filler}.{ext}.log"
#     wildcard_constraints:
#         tool="sickle/se",
#         ext="fastq|fastq.gz"
#     conda:
#         "../envs/sickle.yaml"
#     threads:
#         1
#     shell:
#         "sickle se -f {input.fastq} {params} -o {output.fastq} &> {log}"

