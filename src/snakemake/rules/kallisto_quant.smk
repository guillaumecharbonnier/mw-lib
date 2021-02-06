rule kallisto_quant_pe:
    """
    Created:
        2020-01-30 15:28:50
    Aim:
        Build index for Kallisto.
    Doc:
        https://pachterlab.github.io/kallisto/manual
    Test:
        out/kallisto/quant_pe_kallisto-idx-rnaspades-pe-test-h2al2/sickle/pe_-t_sanger_-q_20/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-KO-Rep1/done
        out/kallisto/quant_pe_kallisto-idx-rnaspades-pe-prc-h2al2/sickle/pe_-t_sanger_-q_20/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-KO-Rep1/abundance.tsv

    out/sickle/pe_-t_sanger_-q_20/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-KO-Rep1_1.fastq.gz
    """
    input:
        fq1="out/{filler}_1.fastq.gz",
        fq2="out/{filler}_2.fastq.gz",
        idx=lambda wildcards: eval(mwconf['ids'][wildcards.kallisto_idx_id])
    output:
        h5="out/{tool}{extra}_{kallisto_idx_id}/{filler}/abundance.h5",
        tsv="out/{tool}{extra}_{kallisto_idx_id}/{filler}/abundance.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="kallisto/quant_pe"
    conda:
        "../envs/kallisto.yaml"
    threads:
        16 #Asking this to try to debug memory shortage.
    shell:
        """
        kallisto quant -i {input.idx} -o `dirname {output.tsv}` {params.extra} {input.fq1} {input.fq2}
        """

