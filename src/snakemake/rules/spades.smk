rule spades_se:
    """
    Created:
        2020-01-20 09:52:58
    Aim:
        Try an alternative to Edena when it fails to run because the input file is to big, e.g. T11C_all_DNA_samples without filtering.
    Test:
        out/spades/se/ln/alias/sst/all_samples/fastq/T11C_H3K27ac/done
    """
    input:
        fq="out/{filler}.fastq.gz"
    output:
        expand("out/{{tool}}{{extra}}/{{filler}}/{filename}",filename=["contigs.fasta", "scaffolds.fasta"]) #Many more files here
    params:
        outdir="out/{tool}{extra}/{filler}",
        extra = params_extra
    wildcard_constraints:
        tool="spades/se"
    conda:
        "../envs/spades.yaml"
    threads:
        16
    shell:
        """
        spades.py -t {threads} {params.extra} -s {input.fq} -o {params.outdir}
        """
