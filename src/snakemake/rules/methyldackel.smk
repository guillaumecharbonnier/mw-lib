rule methyldackel_extract_fagenome:
    """
    Test:
        out/methyldackel/extract_fa-genome-hg19/nextflow/nfcore_methylseq/bwameth/deduplicated/ODG_080.markdup.sorted_CHH.bedGraph
    """
    input:
        bam="out/{filler}.bam",
        bai="out/{filler}.bam.bai",
        fa= lambda wildcards: eval(mwconf['ids'][wildcards.fa_genome_id])
    output:
        bedGraph = "out/{tool}{extra}_{fa_genome_id}/{filler}_{context}.bedGraph"
    params:
        bam      = "out/{tool}{extra}_{fa_genome_id}/{filler}.bam",
        bai      = "out/{tool}{extra}_{fa_genome_id}/{filler}.bam.bai",
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
        ln -srf {input.bam} {params.bam}
        ln -srf {input.bai} {params.bai}
        MethylDackel extract {params.extra} {params.context} {input.fa} {params.bam}
        """


rule tests_methyldackel:
    input:
        expand(
            "out/methyldackel/extract_fa-genome-hg19/samtools/index/nextflow/nfcore_methylseq/bwameth/deduplicated/ODG_{patient}.markdup.sorted_{context}.bedGraph",
            context = ["CpG", "CHG", "CHH"],
            patient = ["080", "081", "082", "083", "084", "085", "086", "087", "088", "089", "090", "091", "092", "093", "094", "095", "096", "097", "098", "099", "100", "101", "102", "103", "104", "105", "106", "107", "108", "109"]
        )



