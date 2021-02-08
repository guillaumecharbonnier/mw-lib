rule bedtools_slop_extra:
    """
    Modified:
        2017-09-29 11:01:06 - add input_chrominfo as chromInfo.
    Doc:
        https://bedtools.readthedocs.io/en/latest/content/tools/slop.html
    Deprecated tests:
       out/bedtools/slop_g-mm10_b-5000/gtftk/5p_3p_coord/sed/add_chr/gunzip/to-stdout/wget/ftp/ftp.ensembl.org/pub/release-87/gtf/mus_musculus/Mus_musculus.GRCm38.87.bed
       out/bedtools/slop_g-hg38_b-1000/gtftk/5p_3p_coord/sed/add_chr/gtftk/select_by_key_key-transcript_id_file-with-values-human-hkg/gunzip/to-stdout/wget/ftp/ftp.ensembl.org/pub/release-91/gtf/homo_sapiens/Homo_sapiens.GRCh38.91.bed
    Test:
        out/bedtools/slop_-b_1000_chrominfo-
    """
    input:
        bed_gff_vcf = "out/{filler}",
        chromInfo = lambda wildcards: eval(mwconf['ids'][wildcards.chrominfo_id])
    output:
        bed_gff_vcf = "out/{tool}{extra}_{chrominfo_id}/{filler}"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="bedtools/slop"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools slop -i {input.bed_gff_vcf} -g {input.chromInfo} {params} > {output.bed_gff_vcf}"

