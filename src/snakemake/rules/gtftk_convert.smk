rule gtftk_convert:
    """
    Created: 
        2017-09-20 12:00:36
    Aim:
        Select only some genes, e.g. coding genes because I do not want to deal with pseudogenes in my heatmaps.
    Doc:
     [-i GTF] [-o GTF] [-n NAME] [-s SEP] [-m MORE_NAMES] [-f {bed,bed3,bed6,bed12}] [-h] [-V ] [-D] [-C] [-K] [-L]

       Description: 
            *  Convert a GTF to various format.

              Version:  2017-09-01

              Arguments:
    -i, --inputfile    Path to the GTF file. Default to STDIN. (default: <stdin>)
    -o, --outputfile   Output file. (default: <stdout>)
    -n, --names        The key(s) that should be used as name. (default: gene_id,transcript_id)
    -s, --separator    The separator to be used for separating name elements (see -n). (default: |)
    -m, --more-names   Add this information to the 'name' column of the BED file. (default: )
    -f, --format       Currently one of bed3, bed6 (default: bed6)

    Test:
        out/gtftk/convert/gtftk/select_by_key_key-gene_id_file-with-values-upreg-ko-nut-r-threshold-1/gunzip/to-stdout/wget/ftp_ensembl/pub/release-87/gtf/mus_musculus/Mus_musculus.GRCm38.87.bed
    """
    input:
        gtf="out/{filler}.gtf"
    output:
        bed="out/gtftk/convert/{filler}.bed"
    conda:
        "../envs/pygtftk.yaml"
    shell:
        "gtftk convert --inputfile {input.gtf} --outputfile {output.bed}"

