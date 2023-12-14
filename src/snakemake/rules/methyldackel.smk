rule methyldackel_extract_fagenome:
    """
    Test:
        out/methyldackel/extract_fa-genome-hg19/samtools/index/nextflow/nfcore_methylseq/bwameth/deduplicated/ODG_081.markdup.sorted_CpG.bedGraph

    -o Output filename prefix. CpG/CHG/CHH context metrics will be output to STR_CpG.bedGraph and so on.
    """
    input:
        bam="out/{filler}.bam",
        bai="out/{filler}.bam.bai",
        fa= lambda wildcards: eval(mwconf['ids'][wildcards.fa_genome_id])
    output:
        bedGraph = "out/{tool}{extra}_{fa_genome_id}/{filler}_{context}.bedGraph"
    params:
        extra = params_extra,
        # function should return "" if context=CpG, "--CHH" if context=CHH and "--CHG" if context=CHG 
        context = lambda wildcards: "--" + wildcards.context if wildcards.context != "CpG" else ""
    wildcard_constraints:
        tool="methyldackel/extract",
        context = "CpG|CHG|CHH"
    conda:
        "../../../../mw-lib/src/snakemake/envs/methyldackel.yaml"
    shell:
        """
        MethylDackel extract {params.extra} {params.context} \
            -o out/{wildcards.tool}{wildcards.extra}_{wildcards.fa_genome_id}/{wildcards.filler} \
            {input.fa} {input.bam}
        """



