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

