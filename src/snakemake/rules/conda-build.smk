rule condabuild:
    """
    Created:
        2018-05-04 11:42:56
    Aim:
        General purpose build rule
    Test:
        out/conda-build/guillaumecharbonnier/chromdet.done
    """
    input:
        conda="opt/miniconda/bin/conda",
        condabuild="opt/miniconda/bin/conda-build",
        meta  = "src/conda/{filler}/meta.yaml",
        build = "src/conda/{filler}/build.sh"
    output:
        done="out/conda-build/{filler}.done"
    shell:
        """
        DIRNAME=`dirname {input.meta}`
        {input.condabuild} $DIRNAME
        touch {output.done}
        """


rule conda_build_gtftk:
    """
    Created:
        2017-02-21 20:52:32
    Modified:
        2017-09-27 15:05:19 - Use directly the conda build and meta from gtftk project.
    Aim:
        Local build of gtftk because I can not upload it on conda cloud.
    """
    input:
        conda="opt/miniconda/bin/conda",
        condabuild="opt/miniconda/bin/conda-build",
        #meta="src/conda/guillaumecharbonnier/gtftk/meta.yaml",
        #build="src/conda/guillaumecharbonnier/gtftk/build.sh"
    output:
        done="src/conda/guillaumecharbonnier/gtftk/done"
    shell:
        """
        # Channel needed for cloud module.
        #{input.conda} config --prepend channels prometeia
        #cd src/conda/guillaumecharbonnier/gtftk
        cd opt/gtftk/conda
        {WDIR}/{input.condabuild} .
        touch {WDIR}/{output.done}
        """


