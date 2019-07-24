rule gtftk_select_by_key_key_value:
    """
    Created: 
        2017-09-20 12:00:36
    Aim:
        Select only some genes, e.g. coding genes because I do not want to deal with pseudogenes in my heatmaps.
    Test:
        out/gtftk/select_by_key_key-gene_biotype_value-protein_coding/gunzip/to-stdout/wget/ftp_ensembl/pub/release-87/gtf/mus_musculus/Mus_musculus.GRCm38.87.gtf
        out/sort/coordinates_bed/gtftk/5p_3p_coord/sed/add_chr/gtftk/select_by_key_key-transcript_biotype_value-protein_coding/gunzip/to-stdout/wget/ftp/ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.bed
    """
    input:
        gtf="out/{filler}.gtf"
    output:
        gtf="out/gtftk/select_by_key_key-{key}_value-{value}/{filler}.gtf"
    wildcard_constraints:
        key="[a-zA-Z0-9-_]+",
        value="[a-zA-Z0-9-_]+"
    conda: 
        "../envs/pygtftk.yaml"
    shell:
        """
        gtftk select_by_key \
            --inputfile {input.gtf} \
            --outputfile {output.gtf} \
            --key {wildcards.key} \
            --value {wildcards.value}
        """

rule gtftk_select_by_key_key_file_with_values:
    """
    Created: 
        2017-09-20 12:00:36
    Aim:
        Select only some genes, e.g. coding genes because I do not want to deal with pseudogenes in my heatmaps.
    Test:
        out/gtftk/select_by_key_key-gene_id_file-with-values-upreg-ko-nut-r-threshold-1/gunzip/to-stdout/wget/ftp_ensembl/pub/release-87/gtf/mus_musculus/Mus_musculus.GRCm38.87.gtf
        input:
            out/sed/gene_name_to_uppercase_in_gtf/gunzip/to-stdout/wget/ftp_ensembl/pub/release-87/gtf/mus_musculus/Mus_musculus.GRCm38.87.gtf
        output:
            expand("out/gtftk/select_by_key_key-gene_name_file-with-values-microarray-ko-nut-threshold-1.5-{pos_neg}/sed/gene_name_to_uppercase_in_gtf/gunzip/to-stdout/wget/ftp_ensembl/pub/release-87/gtf/mus_musculus/Mus_musculus.GRCm38.87.gtf", pos_neg=["pos","neg"])
            out/gtftk/select_by_key_key-gene_name_file-with-values-microarray-ko-nut-threshold-1.5-neg/sed/gene_name_to_uppercase_in_gtf/gunzip/to-stdout/wget/ftp_ensembl/pub/release-87/gtf/mus_musculus/Mus_musculus.GRCm38.87.gtf

            out/gtftk/select_by_key_key-gene_id_file-with-values-gene-id-microarray/gunzip/to-stdout/wget/ftp_ensembl/pub/release-90/gtf/mus_musculus/Mus_musculus.GRCm38.90.gtf

        out/gtftk/select_by_key_key-transcript_id_file-with-values-human-hkg/tar/xvzf_igenome/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf
        out/gtftk/select_by_key_key-transcript_id_file-with-values-human-hkg/gunzip/to-stdout/wget/ftp/ftp.ensembl.org/pub/release-91/gtf/homo_sapiens/Homo_sapiens.GRCh38.91.gtf
        out/gtftk/select_by_key_key-transcript_id_file-with-values-human-hkg/gunzip/to-stdout/wget/ftp/ftp.ensembl.org/pub/release-91/gtf/homo_sapiens/Homo_sapiens.GRCh38.91.chr.gtf
        out/sort/coordinates_bed/gtftk/get_5p_3p_coords/sed/add_chr/gtftk/select_by_key_key-transcript_biotype_file-with-values-txt-gtf-values-protein-coding-TR-IG/gunzip/to-stdout/wget/ftp/ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.bed

    """
    input:
        gtf = "out/{filler}.gtf",
        txt = lambda wildcards: eval(config['ids'][wildcards.txt_id])
    output:
        gtf="out/gtftk/select_by_key_key-{key}_file-with-values-{txt_id}/{filler}.gtf"
    wildcard_constraints:
        key="[a-zA-Z0-9-_]+",
        txt_id="txt-gtf-values-[a-zA-Z0-9-.]+"
    conda:
        "../envs/pygtftk.yaml"
    shell:
        """
        gtftk select_by_key \
            --inputfile {input.gtf} \
            --outputfile {output.gtf} \
            --key {wildcards.key} \
            --file-with-values {input.txt} \
            #--log -V 1 #bugged right now.
        """

