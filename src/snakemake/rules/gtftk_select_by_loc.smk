rule gtftk_select_by_loc_f_n:
    """
    Created: 
        2017-09-20 12:00:36
    Note:
        out/wget/mitra_stanford_ed/kundaje/akundaje/release/blacklists/mm10-mouse/mm10-blacklist.bed.gz
    Aim:
        Select only some genes, e.g. coding genes because I do not want to deal with pseudogenes in my heatmaps.
    Test:
        out/gtftk/select_by_loc_f-mm10-blacklist_invert-match/gunzip/to-stdout/wget/ftp_ensembl/pub/release-87/gtf/mus_musculus/Mus_musculus.GRCm38.87.gtf
    """
    input:
        gtf = "out/{filler}.gtf",
        bed = lambda wildcards: config['ids'][wildcards.bed_id]
    output:
        gtf = "out/gtftk/select_by_loc_f-{bed_id}_invert-match/{filler}.gtf"
    wildcard_constraints:
        bed_id = "[a-zA-Z0-9-]+"
    conda:
        "../envs/pygtftk.yaml"
    shell:
        """
        gtftk select_by_loc \
            --inputfile {input.gtf} \
            --outputfile {output.gtf} \
            --location-file {input.bed} \
            --invert-match
        """

rule gtftk_select_by_loc_f:
    """
    Created: 
        2017-09-20 12:00:36
    Note:
        out/wget/mitra_stanford_ed/kundaje/akundaje/release/blacklists/mm10-mouse/mm10-blacklist.bed.gz
    Aim:
        Select only some genes, e.g. coding genes because I do not want to deal with pseudogenes in my heatmaps.
    Test:
        out/gtftk/select_by_loc_f-mm10-blacklist_invert-match/gunzip/to-stdout/wget/ftp_ensembl/pub/release-87/gtf/mus_musculus/Mus_musculus.GRCm38.87.gtf
    """
    input:
        gtf = "out/{filler}.gtf",
        bed = lambda wildcards: config['ids'][wildcards.bed_id]
    output:
        gtf = "out/gtftk/select_by_loc_f-{bed_id}/{filler}.gtf"
    wildcard_constraints:
        bed_id = "[a-zA-Z0-9-]+"
    conda:
        "../envs/pygtftk.yaml"
    shell:
        """
        gtftk select_by_loc \
            --inputfile {input.gtf} \
            --outputfile {output.gtf} \
            --location-file {input.bed}
        """

