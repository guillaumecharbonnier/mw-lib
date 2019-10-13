rule gtftk_select_most_5p_tx:
    """
    Created: 
        2017-03-09 14:18:23
    Aim:
        Select only the most 5' transcript for each gene, mainly because there are so many very similar transcripts in GRCm38. 
    Doc:
        https://pygtftk.readthedocs.io/en/latest/selection.html#select-most-5p-tx
    Test:
        out/gtftk/select_most_5p_tx/sed/add_chr/gunzip/to-stdout/wget/ftp_ensembl/pub/release-87/gtf/mus_musculus/Mus_musculus.GRCm38.87.gtf
        out/gtftk/select_most_5p_tx/awk/extract_main_chr/sed/add_chr/gunzip/to-stdout/wget/ftp/ftp.ensembl.org/pub/release-95/gtf/mus_musculus/Mus_musculus.GRCm38.95.gtf
    """
    input:
        gtf="out/{filler}.gtf"
    output:
        gtf="out/gtftk/select_most_5p_tx/{filler}.gtf"
    log:
        "out/gtftk/select_most_5p_tx/{filler}.log"
    benchmark:
        "out/gtftk/select_most_5p_tx/{filler}.benchmark.tsv"
    conda:
        "../envs/pygtftk.yaml"
    shell:
        "gtftk select_most_5p_tx --inputfile {input.gtf} > {output.gtf} 2> {log}"

