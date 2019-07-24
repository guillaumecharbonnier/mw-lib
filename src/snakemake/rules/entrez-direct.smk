rule entrez_direct_esearch:
    """
    Test:
        out/entrez-direct/esearch_-db_pubmed_-query_RNA-Seq.txt
    """
    output:
        "out/{tool}{extra}.txt"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "entrez-direct/esearch"
    conda:
        "../envs/entrez-direct.yaml"
    shell:
        "esearch {params.extra} > {output}"
