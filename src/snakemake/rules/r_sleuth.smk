rule r_sleuth_test_h2al2:
    """
    Created:
        2020-01-30 17:04:21
    """
    input:
        expand("out/kallisto/quant_pe_kallisto-idx-rnaspades-pe-test-h2al2/sickle/pe_-t_sanger_-q_20/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-{stage}-H2AL2-{condition}-Rep{replicate}/abundance.tsv", stage=["P","R","C"], condition=["WT","KO"], replicate=["1","2","3"])
    conda:
        "../envs/r_sleuth.yaml"
    shell:
        "echo 'TODO: WRITE SLEUTH RSCRIPT'"
