rule gtftk_retrieve:
    """
    Created: 
        2017-03-07 17:39:30
    Doc:
        https://pygtftk.readthedocs.io/en/latest/information.html#retrieve
    Aim:
        Retrieve gtf file from Ensembl
    Note:
        Unfinished
    Test:
        out/gtftk/retrieve/mus_musculus.gtf.gz
    """
    output:
        gtf="out/{tool}{extra}/{specie}.gtf"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="gtftk/retrieve",
        specie="[a-z]+_[a-z]+"
        #specie="mus_musculus|homo_sapiens|drosophila_melanogaster|saccharomyces_cerevisiae|rattus_norvegicus|[a-zA-Z0-9-]+"
        #specie="^[a-zA-Z][a-z]+_[a-zA-Z][a-z]+|[a-zA-Z0-9-]+" 
        #specie="[mM]us_[mM]usculus|homo_sapiens|[a-zA-Z0-9-]+"
        #specie="mus_musculus|homo_sapiens|[a-z]+_[a-z]+"
    conda:
        "../envs/pygtftk.yaml"
    shell:
        "gtftk retrieve {params.extra} --species-name {wildcards.specie} --to-stdout --delete > {output.gtf}"
        #"gtftk retrieve --species-name {wildcards.species} --to-stdout --delete > {output.gtf}"
        #"gtftk retrieve {params.extra} -o {output.gtf}"
