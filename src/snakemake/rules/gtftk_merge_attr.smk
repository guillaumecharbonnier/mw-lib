rule gtftk_merge_attr:
    """
    Created:
        2017-03-28 18:12:34
    Doc:
        https://pygtftk.readthedocs.io/en/latest/editing.html#merge-attr
    Test:
        out/gtftk/merge_attr_gene_id_gene_name/gunzip/wget/ftp_ensembl/pub/release-87/gtf/mus_musculus/Mus_musculus.GRCm38.87.gtf

        out/gtftk/merge_attr_gene_id_gene_name/cuffmerge/GRCm38_H2AL2_outFilterMultimapNmax-1000/merged.gtf
    """
    input:
        gtf="out/{filler}.gtf"
    output:
        gtf="out/gtftk/merge_attr_gene_id_gene_name/{filler}.gtf"
    conda:
        "../envs/pygtftk.yaml"
    envmodules:
        "pygtftk/1.5.3"
    shell:
        """
        gtftk merge_attr \
            --inputfile {input.gtf} \
            --outputfile {output.gtf} \
            --src-key gene_id,gene_name \
            --dest-key merge_gene_id_name
        """
