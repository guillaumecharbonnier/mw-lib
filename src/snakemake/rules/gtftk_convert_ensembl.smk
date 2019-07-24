rule gtftk_convert_ensembl:
    """
    Test:
        out/gtftk/convert_ensembl/awk/tfbsConsSites_to_gtf.gtf
    """
    input:
        gtf="out/{filler}.gtf"
    output:
        bed="out/gtftk/convert_ensembl/{filler}.gtf"
    conda:
        "../envs/pygtftk.yaml"
    shell:
        "gtftk convert_ensembl --inputfile {input.gtf} --outputfile {output.bed}"
