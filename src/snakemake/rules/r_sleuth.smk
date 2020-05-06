rule r_sleuth_test_h2al2:
    """
    Created:
        2020-01-30 17:04:21
    Note:
        Extract transcripts containing GSAT using:
        F=out/repeatmasker/test2_-species_mouse_-a_-html/sed/trim-fasta-header-from-rnaspades/gunzip/to-stdout/ln/updir/mw-sk/inp/transcripts/transcripts.fasta.out;  head -2 $F; grep -e "NODE.*GSAT" $F | awk '{print $5}'
    
    out/kallisto/quant_pe_kallisto-idx-rnaspades-pe-prc-h2al2/sickle/pe_-t_sanger_-q_20/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-KO-Rep1/abundance.tsv
    Note:
        here input gather multiple transcriptome assemblies with differentes names for transcript NODE[0-9]+
    """
    input:
        quant=expand("out/kallisto/quant_pe_-b_10_kallisto-idx-rnaspades-pe-prc-h2al2/sickle/pe_-t_sanger_-q_20/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-{stage}-H2AL2-{condition}-Rep{replicate}/abundance.tsv", stage=["P","R","C"], condition=["WT","KO"], replicate=["1","2","3"], index=["kallisto-idx-rnaspades-pe-test5-h2al2","kallisto-idx-rnaspades-pe-test-h2al2","kallisto-idx-rnaspades-pe-test4-h2al2","kallisto-idx-rnaspades-pe-test3-c-h2al2","kallisto-idx-rnaspades-pe-test2-c-h2al2-ko","kallisto-idx-rnaspades-pe-test1-h2al2"]),
        rmsk="out/repeatmasker/test2_-species_mouse_-a_-html/sed/trim-fasta-header-from-rnaspades/gunzip/to-stdout/ln/updir/mw-sk/inp/transcripts/transcripts.fasta.out"
    conda:
        "../envs/r_sleuth.yaml"
    shell:
        "echo 'TODO: WRITE SLEUTH RSCRIPT'"
