rule fastx_trimmer_legacy:
    """
    Created:
        2016-07-13 11h38
    Modified:
        2016-11-22 15h15 - New filler pattern.
        2017-05-06 14:16:55 - Updated to conda package.
    Aim:
        Trim fastq to a fixed length. Used before calculation for insert size to avoid artefacts.
    Note:
        You need to add the -Q33 parameter to tell it that you're using Sanger encoded quality scores, not obsolete Illumina one.
    Test:
        out/fastx_toolkit/fastx_trimmer_l-30/sra-tools/fastq-dump_se/SRR2096391.fastq
        out/fastx_toolkit/fastx_trimmer_l-30/sra-tools/fastq-dump_se/SRR646457.fastq.gz
        
    """   
    input:
        "out/{filler}.fast{aq}{gz}",
    output:
        "out/fastx_toolkit/fastx_trimmer_l-{l}/{filler}.fast{aq}{gz}"
    params:
        gz = lambda wildcards: "-z" if wildcards.gz == '.gz' else "",
        cat = lambda wildcards: "zcat" if wildcards.gz == '.gz' else "cat"
    log:
        "out/fastx_toolkit/fastx_trimmer_l-{l}/{filler}.fast{aq}{gz}.log"
    benchmark:
        "out/fastx_toolkit/fastx_trimmer_l-{l}/{filler}.fast{aq}{gz}.benchmark.tsv"
    wildcard_constraints:
        l="[0-9]+",
        aq="a|q",
        gz=".gz|"
    conda:
        "../envs/fastx_toolkit.yaml"
    shell:
        "{params.cat} {input} | fastx_trimmer -l {wildcards.l} -Q33 {params.gz} -o {output} &> {log}"

