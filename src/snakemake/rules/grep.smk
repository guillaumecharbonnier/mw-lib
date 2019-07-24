rule grep_extra:
    """
    Created:
        2019-02-01 11:50:39
    Aim:
    Test:
        out/grep/ENSG00000277734/sed/add_chr/gunzip/to-stdout/wget/ftp/ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.gtf
    """
    input:
        txt="out/{filler}"
    output:
        txt="out/{tool}{extra}/{filler}"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="grep/"
    conda:
        "../envs/coreutils.yaml"
    shell:
        "grep {params.extra} {input.txt} > {output.txt}"
