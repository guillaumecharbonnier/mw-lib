rule kallisto_index:
    """
    Created:
        2020-01-30 15:28:50
    Aim:
        Build index for Kallisto.
    Doc:
        https://pachterlab.github.io/kallisto/manual
    Test:
        out/kallisto/index/rnaspades/pe_test_h2al2/transcripts.idx

    """
    input:
        "out/{filler}.fasta"
    output:
        "out/{tool}{extra}/{filler}.idx"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="kallisto/index"
    conda:
        "../envs/kallisto.yaml"
    shell:
        """
        kallisto index {params.extra} -i {output} {input}
        """
