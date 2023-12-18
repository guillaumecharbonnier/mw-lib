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
        tool = "bismark/methylation_extractor",
    shell:
        """
        # Seems we need to move to output directory else bismark_methylation_extractor will flood the working directory with output files
        # Moreover, since the same filename can come from different directories, we need to make sure we don't overwrite files
        WDIR=`pwd`
        # Move to the output directory
        cd `dirname {output.done}`
        bismark_methylation_extractor {params.extra} $WDIR/{input.bam} &> {log}
        touch $WDIR/{output.done}
        """
