rule pbgzip_tabix_bed:
    """
    Created:
        2019-10-15 18:14:50
    Aim:
        Test the combination of pbgzip and tabix to compress and index bed files.
        Originally done to test if Cyverse can recognize automatically such bed compressed files as bedgz.
    Test:

    """
    input:
        bed="out/{filler}.bed",
    output:
        bedgz="out/pbgzip_tabix_bed/{filler}.bed.gz",
        tbi="out/pbgzip_tabix_bed/{filler}.bed.gz.tbi"
    params:
        bed="out/pbgzip_tabix_bed/{filler}.bed"
    log:
            "out/pbgzip_tabix_bed/{filler}.log"
    benchmark:
            "out/pbgzip_tabix_bed/{filler}.benchmark.tsv"
    conda:
        "../envs/pbgzip_tabix.yaml"
    threads:
        1
    shell:
        """
        (ln -srf {input.bed} {params.bed}
        pbgzip {params.bed}
        tabix -p bed {output.bedgz}
        ) &> {log}
        """

