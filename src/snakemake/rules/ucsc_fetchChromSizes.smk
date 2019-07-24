rule ucsc_fetchChromSizes:
    """
    Created:
        2018-03-06 18:31:08
    Doc:
        http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes
    Test:
        out/ucsc/fetchChromSizes/hg38.txt
    """
    output:
        "out/ucsc/fetchChromSizes/{assembly}.txt"
    conda:
        "../envs/ucsc.yaml"
    shell:
        "fetchChromSizes {wildcards.assembly} > {output}"
