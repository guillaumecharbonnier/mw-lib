rule r_tidy_featureCounts:
    """
    Aim:
        Tidy featureCounts output to a simple matrix.
        raw.tsv can notably be used as input to:
        https://amp.pharm.mssm.edu/biojupies/
    Test:
        out/r/tidy_featureCounts/subread/featureCounts_-O_gtf-GRCh38-ensembl_bam-GRCh38-tgml-rna-kit-benchmark_raw.tsv
    """
    input:
        tsv = "out/{filler}.tsv"
    output:
        raw = "out/r/tidy_featureCounts/{filler}_raw.tsv",
        rpkm = "out/r/tidy_featureCounts/{filler}_rpkm.tsv"
    conda:
        "../envs/r_tidy_featureCounts.yaml"
    script:
        "../../r/scripts/tidy_featureCounts.R"
