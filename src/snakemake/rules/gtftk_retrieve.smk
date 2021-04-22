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
        gtf="out/{tool}{extra}/{specie}.gtf.gz"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="gtftk/retrieve",
        specie="[mM]us_[mM]usculus|homo_sapiens|[dD]rosophila_[mM]elanogaster|[sS]accharomyces_[cC]erevisiae|[rR]attus_[nN]orvegicus|[a-zA-Z0-9-]+"
        #specie="^[a-zA-Z][a-z]+_[a-zA-Z][a-z]+|[a-zA-Z0-9-]+" 
        #specie="[mM]us_[mM]usculus|homo_sapiens|[a-zA-Z0-9-]+"
        #specie="mus_musculus|homo_sapiens|[a-z]+_[a-z]+"
    conda:
        "../envs/pygtftk.yaml"
    shell:
        "gtftk retrieve {params.extra} -o {output.gtf}"
        #"gtftk retrieve --species-name {wildcards.species} --to-stdout --delete > {output.gtf}"

