rule bbduk_se:
    """
    Created:
    2020-08-18 15:06:17
    Usage:
        bbduk.sh in=<input> out=<reads to keep> outt=<reads to toss> hist=<histogram output>
    IUse in2 for paired reads in a second file  
    Test:
        out/bbduk/se/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac.fastq.gz
    """
    input:
        fq="out/{filler}.fastq.gz",
    output:
        fq_to_keep="out/{tool}{extra}/{filler}.fastq.gz"#,
        #fq_to_toss="out/{tool}{extra}/{filler}.tossed.fastq.gz"
    log:
        "out/{tool}{extra}/{filler}.log"
    benchmark:
        "out/{tool}{extra}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "bbduk/se",
    conda:
        "../envs/bbmap.yaml"
    shell:
        "bbduk.sh {params.extra} in='{input.fq}' out='{output.fq_to_keep}' &> {log}"

rule bbduk_pe:
    """
    Created:
        2020-08-18 15:06:17
    Usage:
        bbduk.sh in=<input> out=<reads to keep> outt=<reads to toss> hist=<histogram output>
    IUse in2 for paired reads in a second file  
    Test:
        out/bbduk/pe/ln/alias/sst/all_samples/fastq/856_H3K27ac_1.fastq.gz
    """
    input:
        fq_1="out/{filler}_1.fastq.gz",
        fq_2="out/{filler}_2.fastq.gz"
    output:
        fq_1 = "out/{tool}{extra}/{filler}_1.fastq.gz",
        fq_2 = "out/{tool}{extra}/{filler}_2.fastq.gz"
        #fq_to_toss="out/{tool}{extra}/{filler}.tossed.fastq.gz"
    log:
        "out/{tool}{extra}/{filler}.log"
    benchmark:
        "out/{tool}{extra}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "bbduk/pe"
    conda:
        "../envs/bbmap.yaml"
    shell:
       "bbduk.sh {params.extra} in='{input.fq_1}' in2='{input.fq_2}' out='{output.fq_1}' out2='{output.fq_2}' &> {log}"



