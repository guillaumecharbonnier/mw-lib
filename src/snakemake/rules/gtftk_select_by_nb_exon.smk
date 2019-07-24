rule gtftk_select_by_nb_exon_m_M:
    """
    Created:
        2017-09-19 17:17:07
    Aim:
        Select transcripts based on the number of exons.
    Doc:
        https://pygtftk.readthedocs.io/en/latest/selection.html#select-by-nb-exon
    Test:
        out/gtftk/select_by_nb_exon_m-2_M-3/{filler}.gtf
    """
    input:
        gtf="out/{filler}.gtf"
    output:
        gtf="out/gtftk/select_by_nb_exon_m-{min}_M-{max}/{filler}.gtf"
    conda:
        "../envs/pygtftk.yaml"
    shell:
        """
        gtftk select_by_nb_exon \
            --inputfile {input.gtf} \
            --outputfile {output.gtf} \
            -m {output.min} \
            -M {output.max}
        """
        

