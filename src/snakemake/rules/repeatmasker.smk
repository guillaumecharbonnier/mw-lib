rule repeatmasker_test1:
    """
    out/repeatmasker/test1_-species_mouse_-a_-html/rnaspades/pe_test_h2al2/transcripts/done

    """
    input:
        fasta="out/{filler}.fasta"
    output:
        touch("out/{tool}{extra}/{filler}/done")
    params:
        extra = params_extra
    wildcard_constraints:
        tool="repeatmasker/test1"
    conda:
        "../envs/repeatmasker.yaml"
    threads:
        16
    shell:
        """
        RepeatMasker -parallel {threads} {params.extra} --dir `dirname {output}` {input.fasta}
        """
