rule gtftk_short_long:
    """
    Created:
        2017-09-19 16:23:38
    Aim:
        By default, short_long will output only shortest transcript for each gene.
    Doc:
        https://pygtftk.readthedocs.io/en/latest/selection.html#short-long
    Test:
        out/gtftk/short_long/gtftk/merge_attr_gene_id_gene_name/gunzip/to-stdout/wget/ftp_ensembl/pub/release-87/gtf/mus_musculus/Mus_musculus.GRCm38.87.gtf
    """
    input:
        gtf="out/{filler}.gtf"
    output:
        gtf="out/gtftk/short_long/{filler}.gtf"
    conda:
        "../envs/pygtftk.yaml"
    shell:
        "gtftk short_long --inputfile {input.gtf} --outputfile {output.gtf}"
        
rule gtftk_short_long_longs:
    """
    Created:
        2017-09-19 16:29:02
    Aim:
        With --longs, short_long will output only longest transcript for each gene.
    Test:
        out/gtftk/short_long_longs/gtftk/merge_attr_gene_id_gene_name/gunzip/to-stdout/wget/ftp_ensembl/pub/release-87/gtf/mus_musculus/Mus_musculus.GRCm38.87.gtf
    """
    input:
        gtf="out/{filler}.gtf"
    output:
        gtf="out/gtftk/short_long_longs/{filler}.gtf"
    conda:
        "../envs/pygtftk.yaml"
    shell:
        "gtftk short_long --longs --inputfile {input.gtf} --outputfile {output.gtf}"
