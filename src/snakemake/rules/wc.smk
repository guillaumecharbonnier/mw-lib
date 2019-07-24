rule wc_extra:
    """
    Created:
        2017-05-20 10:54:44
    Aim:
    Test:
        out/wc/_-l/sra-tools/fastq-dump_se/SRR1202037.fastq.txt
        out/wc/_/sra-tools/fastq-dump_se/SRR1202037.fastq.txt
        
    """
    input:
        txt="out/{filler}"
    output:
        txt="out/{tool}{extra}/{filler}.txt"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="wc/"
    shell:
        "wc {params.extra} {input.txt} > {output.txt}"

