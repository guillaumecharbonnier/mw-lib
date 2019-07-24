rule cp_alias:
    """
    Aim:
        Alternative to ln_alias. Is useful in case where alias target is write-protected because snakemake will currently refuse to process them. Probably a bug.
    Test:
    """
    input:
        lambda wildcards: config['ids'][wildcards.id]
        #input_ln_alias
    output:
        alias="out/cp/alias/{id}"
    #wildcard_constraints:
    #    alias_id=".*"
    conda:
        "../envs/coreutils.yaml"
    shell:
        # quotes around {input} is mandatory to be able
        # to have files with special characters in input.
        "cp '{input}' {output.alias}"
