rule umi_tools_dedup:
    """
    Created:
        2020-02-03 18:35:16
    Doc:
        https://umi-tools.readthedocs.io/en/latest/reference/dedup.html
    Note:
        Compare with https://raw.githubusercontent.com/e-hutchins/dqRNASeq/master/dqRNASeq.sh
    TODO:
        Compare multiple treatments with umitools to understand better what is happening in TGML report.
    Test:
        out/samtools/index/umi-tools/dedup_--method_unique/star/pe_fastq.gz_to_bam_standard_staridx-GRCh38-ensembl_gtf-GRCh38-ensembl/sickle/pe_-t_sanger_-q_20/umi-tools/extract_pe_--bc-pattern=NNNNNNNN_--bc-pattern2=NNNNNNNN/ln/alias/fastq/UMI_Rg_GFP.bam
    """
    input:
        bam="out/{filler}.bam",
        bai="out/{filler}.bam.bai"
    output:
        bam="out/{tool}{extra}/{filler}.bam",
        tsv=expand("out/{{tool}}{{extra}}/{{filler}}_{tsv_suffix}.tsv", tsv_suffix=["per_umi_per_position","per_umi", "edit_distance"])
    log:
            "out/{tool}{extra}/{filler}.log"
    benchmark:
            "out/{tool}{extra}/{filler}.benchmark.tsv"
    params:
        extra = params_extra,
        tsv_prefix="out/{tool}{extra}/{filler}"
    wildcard_constraints:
        tool="umi-tools/dedup"
    conda:
        "../envs/umi_tools.yaml"
    shell:
        """
        umi_tools dedup -I {input.bam} --output-stats {params.tsv_prefix} -S {output.bam} {params.extra} -L {log}
        """


rule dev_umi_tools_dedup_discrepancy_investigation:
    """
    TODO:
        Produce these files then look in IGV to the first exon of TMSB4X.
    """
    input:
        expand("out/samtools/index/ln/alias/dev_umi_tools_dedup_discrepancy_investigation_Rg_GFP/{sample}.{ext}", sample=["UMI_ND", "Std_ND", "UMI_STL", "UMI_ST2", "UMI_ST3"], ext=["bam", "bam.bai"])

