rule gtftk_get_5p_3p_coords:
    """
    Created: 
        2017-08-03 13:01:22
    Modified:
        2018-10-31 18:26:38 - Function is now called "get_5p_3p_coord" and should produce outputs in a path accordingly to it.
    Aim:
        Select only the most 5' transcript for each gene, mainly because there are so many very similar transcripts in GRCm38. 
    Doc:
        https://pygtftk.readthedocs.io/en/latest/selection.html#short-long
    Test:
        out/gtftk/get_5p_3p_coords/sed/add_chr/gunzip/to-stdout/wget/ftp/ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.bed
        out/gtftk/get_5p_3p_coords_-n_gene_name/sed/add_chr/gunzip/to-stdout/wget/ftp/ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.bed
        out/gtftk/get_5p_3p_coords_-n_gene_name/sed/add_chr/gunzip/to-stdout/wget/ftp/ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.bed
    """
    input:
        gtf="out/{filler}.gtf"
    output:
        bed="out/{tool}{extra}/{filler}.bed"
    wildcard_constraints:
        tool="gtftk/get_5p_3p_coords"
    params:
        extra=params_extra
    conda:
        "../envs/pygtftk.yaml"
    shell:
        "gtftk get_5p_3p_coords {params.extra} --inputfile {input.gtf} --outputfile {output.bed}"

