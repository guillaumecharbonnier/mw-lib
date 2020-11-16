rule tadpole_se:
    """
    Created:
        2020-08-18 15:06:17
    Aim:
        tadpole.sh mode=correct in=reads.fq out=corrected.fq

        tadpole.sh mode=correct in1=read1.fq in2=read2.fq out1=corrected1.fq out2=corrected2.fq

        https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/tadpole-guide/
        supports every value of K from 1-31, every multiple of 2 from 32-62 (meaning 32, 34, 36, etc), every multiple of 3 from 63-93
        Typically, about 2/3rds of read length is a good value for K for assembly.

    Test:
        out/tadpole/assembly_se/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac.fasta
    """
    input:
        fq="out/{filler}.fastq.gz",
    output:
        fa="out/{tool}{extra}/{filler}.fasta"
    log:
        "out/{tool}{extra}/{filler}.log"
    benchmark:
        "out/{tool}{extra}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "tadpole/assembly_se",
    conda:
        "../envs/bbmap.yaml"
    shell:
        "tadpole.sh {params.extra} in='{input.fq}' out='{output.fa}' &> {log}"
