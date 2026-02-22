rule mixcr_analyze_pe:
    """
    Created:
        2025-05-26
    Doc:
        https://mixcr.readthedocs.io/en/master/
    Aim:
        Aligns sequencing reads to V(D)J reference using MiXCR.
    Test:
        out/mixcr/analyze_rna-seq_--species_hsa/agent/trim_-v2/ln/alias/sst/all_samples/fastq/1379_BOUNAD_RNA_XT_HS2.clna
    """
    input:
        fq1 = "out/{filler}_1.fastq.gz",
        fq2 = "out/{filler}_2.fastq.gz"
        # ref = lambda wildcards: eval(mwconf['ids'][wildcards.mixcr_ref_id])
    output:
        # clna = "out/{tool}{extra}_{mixcr_ref_id}/{filler}.clna"
        clna = "out/{tool}{extra}/{filler}.clna"
    # log:
    #     "out/{tool}{extra}/{filler}.log"
        # "out/{tool}{extra}_{mixcr_ref_id}/{filler}.log"
    benchmark:
        # "out/{tool}{extra}_{mixcr_ref_id}/{filler}.benchmark.tsv"
        "out/{tool}{extra}/{filler}.benchmark.tsv"
    params:
        extra = params_extra,
        outprefix = "out/{tool}{extra}/{filler}"
    wildcard_constraints:
        tool = "mixcr/analyze_pe"
    conda:
        "../envs/mixcr.yaml"
    priority: 10
    threads: 2
    shell:
        """
        # mixcr analyze {params.extra} -f {input.fq1} {input.fq2} {params.outprefix}
        # Replacing the above command with the following to try to avoid these kinds of issues:
        # https://github.com/guillaumecharbonnier/mw-lib/issues/38
        tmpdir=$(mktemp -d /tmp/mixcr_analyze_pe.XXXXXX)
        tmpprefix=$(basename {params.outprefix})
        mkdir -p "$tmpdir"/"$tmpprefix"
        mixcr analyze {params.extra} -f {input.fq1} {input.fq2} "$tmpdir"/"$tmpprefix"
        mv "$tmpdir"/* $(dirname {params.outprefix})/
        rmdir "$tmpdir"
        """

rule mixcr_analyze_se:
    """
    Created:
        2025-05-26
    Doc:
        https://mixcr.readthedocs.io/en/master/
    Aim:
        Aligns sequencing reads to V(D)J reference using MiXCR.
    Test:
        Never tried
    """
    input:
        fq = "out/{filler}.fastq.gz",
    output:
        clna = "out/{tool}{extra}/{filler}.clna"
    # log:
    #     "out/{tool}{extra}/{filler}.log"
        # "out/{tool}{extra}_{mixcr_ref_id}/{filler}.log"
    benchmark:
        # "out/{tool}{extra}_{mixcr_ref_id}/{filler}.benchmark.tsv"
        "out/{tool}{extra}/{filler}.benchmark.tsv"
    params:
        extra = params_extra,
        outdir = "out/{tool}{extra}/{filler}"
    wildcard_constraints:
        tool = "mixcr/analyze_se"
    conda:
        "../envs/mixcr.yaml"
    priority: 10
    threads: 2
    shell:
        """
        mixcr analyze {params.extra} -f {input.fq} {params.outdir}
        """