#http://gage.cbcb.umd.edu/results/index.html
#https://github.com/bcgsc/abyss
rule abyss_pe:
    """
    Created:
        2019-12-19 17:50:15
    Test:
        out/abyss/pe_k=30_fastq.gz/ln/alias/sst/all_samples/fastq/489_H3K27ac-contigs.fa out/abyss/pe_k=25_fastq.gz/ln/alias/sst/all_samples/fastq/489_H3K27ac-contigs.fa
    """
    input:
        mate1="out/{filler}_1.{ext}",
        mate2="out/{filler}_2.{ext}"
    output:
        "out/{tool}{extra}_{ext}/{filler}-contigs.fa"
    log:
        "out/{tool}{extra}_{ext}/{filler}/log"
    benchmark:
        "out/{tool}{extra}_{ext}/{filler}/benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "abyss/pe",
        ext="fa|fq|fasta|fastq|fasta.gz|fastq.gz"
    threads:
        16
    conda:
        "../envs/abyss.yaml"
    shell:
        "abyss-pe {params.extra} np={threads} name=out/{wildcards.tool}{wildcards.extra}_{wildcards.ext}/{wildcards.filler} in='{input}' &> {log}"

rule abyss_debug_pe:
    """
    Created:
        2019-12-19 17:50:15
    Test:
        out/abyss/debug_pe_fastq.gz/bbmap/reformat_fastq_pe_addslash/ln/alias/sst/all_samples/fastq/489_H3K27ac-contigs.fa
    """
    input:
        mate1="out/{filler}_1.{ext}",
        mate2="out/{filler}_2.{ext}"
    output:
        "out/{tool}{extra}_{ext}/{filler}-contigs.fa"
    log:
        "out/{tool}{extra}_{ext}/{filler}/log"
    benchmark:
        "out/{tool}{extra}_{ext}/{filler}/benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "abyss/debug_pe",
        ext="fa|fq|fasta|fastq|fasta.gz|fastq.gz"
    threads:
        16
    conda:
        "../envs/abyss.yaml"
    shell:
        "abyss-pe {params.extra} k=25 np={threads} name=out/{wildcards.tool}{wildcards.extra}_{wildcards.ext}/{wildcards.filler} in='{input}' &> {log}"

