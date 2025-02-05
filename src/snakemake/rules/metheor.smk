rule metheor_extra:
    """
    Created:
        2025-01-13 11:00:19
    Doc:
        https://github.com/dohlee/metheor
    Note:
    Test:
        out/metheor/me/fumi_tools/dedup_--paired/samtools/sort/bismark/pe_--pbat_fa-genome-hs1/trim_galore/pe_--non_directional_--rrbs_--length_15_--fastqc/ln/updir/mw-tall-rrbs/inp/bcl/240521_A00680_0489_BHVC2MDMXY/data-isilon/raw-data/novaseq_imagine/240521_A00680_0489_BHVC2MDMXY/Data/Intensities/BaseCalls/fumi_tools_demultiplex/477_DOUYAN.tsv
    """
    input:
        bam="out/{filler}.bam",
        # bai="out/{filler}.bam.bai"#,
        # fa= lambda wildcards: eval(mwconf['ids'][wildcards.fa_genome_id])
    output:
        tsv="out/{tool}{extra}/{filler}.tsv"
    log:
            "out/{tool}{extra}/{filler}.log"
    benchmark:
            "out/{tool}{extra}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    conda:
        "../envs/metheor.yaml"
    wildcard_constraints:
        tool="metheor/"
    # threads:
    #     MAX_THREADS
    shell:
        """
        metheor {params.extra} --input {input.bam}  --output {output.tsv} &> {log}
        """
