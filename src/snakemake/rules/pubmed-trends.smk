rule pubmed_trends:
    """
    Uncompleted
    Test:
        out/academic-keyword-occurence/RNA-Seq_2000_2018.csv
    """
    output:
        "out/{tool}{extra}.csv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "academic-keyword-occurence/"
    conda:
        "../envs/academic-keyword-occurence.yaml"
    shell:
        "extract_occurrences.py {params.extra} {output}"
