rule wiggletools_wigToBigWig_extra:
    """
    Created:
        2013-02-14 
    Aim:
        Operations on all wig/bw/... in a folder, followed by a conversion from wig to bigwig.
    Note:
        -clip is nearly mandatory to solve some issues.
    Test:

    """
    input:
        folder = "out/{filler}",
        chrominfo = lambda wildcards: eval(mwconf['ids'][wildcards.chrominfo_id])
    output:
        bw = "out/{tool}_{chrominfo_id}/{filler}/{operator}.bw"
        # bw = "out/{tool}{extra}{chrominfo_id}/{filler}/{operator}.bw"
    # params:
    #     extra = params_extra
    wildcard_constraints:
        tool="wiggletools_wigToBigWig/from_folder",
        operator="sum|median|mean|min|max" # Other operators could be added
    priority:
        3
    conda:
        "../envs/wiggletools.yaml"
    shell:
        """
        cd {input.folder}
        wiggletools {wildcards.operator} * | wigToBigWig stdin {WDIR}/{input.chrominfo} {WDIR}/{output.bw}
        """
