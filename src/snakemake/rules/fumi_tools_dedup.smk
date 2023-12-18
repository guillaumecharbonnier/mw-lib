rule fumi_tools_dedup:
    """
    Created:
        2023-12-18 18:16:09
    Note:
        Sort bam by read name should be used for example before bedtools bamtobed bedpe.
    Test:

    """
    input:
        bam="out/{filler}.bam"
    output:
        bam="out/{tool}{extra}/{filler}.bam"
    log:
        "out/{tool}{extra}/{filler}.log"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="fumi_tools/dedup"
    threads:
        MAX_THREADS
    # The version on conda is not the latest one
    # Relying on manual system install from 
    # https://ccb-gitlab.cs.uni-saarland.de/tobias/fumi_tools
    # for now
    # conda:
    #     "../envs/samtools.yaml"
    # NOTE THAT fumi_tools requires samtools to be installed for this step
    # and python
    shell:
        """
        fumi_tools dedup -i {input.bam} -o {output.bam} {params.extra} --threads {threads} &> {log}
        """

