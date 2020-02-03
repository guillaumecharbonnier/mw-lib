rule umi_tools_dedup:
    """
    Created:
        2020-02-03 18:35:16
    Note:
        Compare with https://raw.githubusercontent.com/e-hutchins/dqRNASeq/master/dqRNASeq.sh
    Test:
        out/umi-tools/dedup/star/pe_fastq.gz_to_bam_standard_staridx-GRCh38-ensembl_gtf-GRCh38-ensembl/sickle/pe_-t_sanger_-q_20/umi-tools/extract_pe_--bc-pattern=NNNNNNNN_--bc-pattern2=NNNNNNNN/ln/alias/fastq/UMI_RgPML_3K.bam
    """
    input:
        bam="out/{filler}.bam",
        bai="out/{filler}.bam.bai"
    output:
        bam="out/{tool}{extra}/{filler}.bam",
        tsv="out/{tool}{extra}/{filler}.tsv"
    log:
            "out/{tool}{extra}/{filler}.log"
    benchmark:
            "out/{tool}{extra}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="umi-tools/dedup"
    conda:
        "../envs/umi_tools.yaml"
    shell:
        """
        umi_tools dedup -I {input.bam} --output-stats={output.tsv} -S {output.bam} {params.extra} -L {log}
        """


