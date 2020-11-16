rule umi_tools_extract_se:
    """
    Created:
        2020-09-04 16:52:29
    Regex for Nextflex miRNA v3:
        ^(?P<umi_1>.{4}).+(?P<umi_2>.{4})$
        https://github.com/CGATOxford/UMI-tools/issues/356
    Test:
        out/umi-tools/extract_se_NEXTflex-Small-RNA-seq-v3/cutadapt/se_-a_TGGAATTCTCGGGTGCCAAGG_--minimum-length_23/ln/updir/mw-el-cherif/inp/fastq/Run_295/S003257_batch_A_297_603-1_R1.fastq.gz
    """
    input:
        fq="out/{filler}.{ext}"
    output:
        fq="out/{tool}{extra}/{filler}.{ext}"
    log:
           "out/{tool}{extra}/{filler}.{ext}.log"
    benchmark:
           "out/{tool}{extra}/{filler}.{ext}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="umi-tools/extract_se",
        ext="fastq|fq|fastq.gz|fq.gz"
    conda:
        "../envs/umi_tools.yaml"
    shell:
        """
        umi_tools extract -I {input.fq} --stdout={output.fq} {params.extra} -L {log}
        """

rule umi_tools_extract_pe:
    """
    Created:
        2020-02-03 15:57:24
    Test:
        out/umi-tools/extract_pe_--bc-pattern=NNNNNNNN_--bc-pattern2=NNNNNNNN/ln/alias/fastq/UMI_RgPML_3K_1.fastq.gz
    """
    input:
        R1="out/{filler}_1.{ext}",
        R2="out/{filler}_2.{ext}"
    output:
        R1="out/{tool}{extra}/{filler}_1.{ext}",
        R2="out/{tool}{extra}/{filler}_2.{ext}"
    log:
           "out/{tool}{extra}/{filler}.{ext}.log"
    benchmark:
           "out/{tool}{extra}/{filler}.{ext}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="umi-tools/extract_pe",
        ext="fastq|fq|fastq.gz|fq.gz"
    conda:
        "../envs/umi_tools.yaml"
    shell:
        """
        umi_tools extract -I {input.R1} --read2-in={input.R2} --stdout={output.R1} --read2-out={output.R2} {params.extra} -L {log}
        """

