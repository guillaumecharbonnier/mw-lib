rule bismark_methylation_extractor:
    input:
        bam = "out/{filler}.bam"
    output:
        done = "out/{tool}{extra}/{filler}.done"
    log:
        "out/{tool}{extra}/{filler}.log"
    params:
        extra = params_extra
    conda:
        "../envs/bismark.yaml"
    wildcard_constraints:
        tool = "bismark/methylation_extractor"
    shell:
        """
        bismark_methylation_extractor {params.extra} -o `dirname {output.done}` {input.bam} &> {log}
        touch {output.done}
        """
