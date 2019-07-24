localrules: stow_extra
rule stow_extra:
    """
    Created:
        2018-11-07 01:54:53
    Test:
        out/stow/r/associate_atac_peaks_to_nearest_genes/done
    """
    input:
        lambda wildcards: config['ids'][wildcards.id]
    output:
        "out/stow/{id}"
        #"out/{tool}{extra}/{filler}/done"
    #wildcard_constraints:
    #    tool="stow/"
    #params:
    #    extra = params_extra
    resources:
        stow_token = 1
    conda:
        "../envs/stow.yaml"
    shell:
        "INPUTDIR=`dirname {input}` ;"
        "stow -d `dirname $INPUTDIR` -t `dirname {output}` `basename $INPUTDIR`"
