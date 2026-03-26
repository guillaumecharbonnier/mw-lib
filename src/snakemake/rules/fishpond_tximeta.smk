"""
Aim:
    Rules to import salmon quant.sf files with tximeta (human GRCh38 / Ensembl)
    and run fishpond preprocessing, producing:
      - transcript-level SummarizedExperiment (transcripts.rds)
      - gene-level SummarizedExperiment      (genes.rds)

Usage:
    mapper.py populates mwconf['ids'][coldata_id] with a list of
    "sample_name\\tquant_sf_path" strings for each RNA exp group (GRCh38).
    The rules below write the coldata TSV and run the R script.

    Alternatively, register a list of "name\\tfile" strings manually in
    mwconf['ids'] under a chosen <coldata_id> key and request:
        out/fishpond/tximeta_GRCh38_ensembl/<coldata_id>/transcripts.rds
        out/fishpond/tximeta_GRCh38_ensembl/<coldata_id>/genes.rds
"""


rule fishpond_coldata_from_ids:
    """
    Created:
        2026-02-22
    Aim:
        Write a coldata TSV (columns: name, file) from the list of
        "name\\tfile" strings stored in mwconf['ids'][coldata_id].
        Populated automatically by mapper.py for each RNA exp group.
    """
    output:
        "out/fishpond/coldata/{coldata_id}.tsv"
    wildcard_constraints:
        coldata_id = r"[-a-zA-Z0-9_]+"
    run:
        import os
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        with open(output[0], 'w') as fh:
            fh.write("name\tfile\n")
            for entry in mwconf['ids'][wildcards.coldata_id]:
                fh.write(entry + "\n")


rule fishpond_tximeta_GRCh38_ensembl:
    """
    Created:
        2026-02-22
    Aim:
        Import salmon quant.sf files (human GRCh38, Ensembl release detected
        automatically from salmon index metadata) with tximeta, add SYMBOL IDs,
        scale inferential replicates, and summarise to gene level.
    coldata_id:
        Key in mwconf['ids'] whose value is a list of "name\\tfile" strings.
        Populated automatically by mapper.py for RNA exp groups.
    Test:
        out/fishpond/tximeta_GRCh38_ensembl/coldata-GRCh38-myexp/transcripts.rds
        out/fishpond/tximeta_GRCh38_ensembl/coldata-GRCh38-myexp/genes.rds
        out/fishpond/tximeta_GRCh38_ensembl/coldata-GRCh38-myexp/transcripts_counts.csv.gz
        out/fishpond/tximeta_GRCh38_ensembl/coldata-GRCh38-myexp/genes_counts.csv.gz
    """
    input:
        coldata  = "out/fishpond/coldata/{coldata_id}.tsv",
        quant_sf = lambda wildcards: [e.split('\t')[1] for e in mwconf['ids'][wildcards.coldata_id]],
        script   = "../mw-lib/src/r/scripts/fishpond_tximeta.R"
    output:
        transcripts        = "out/fishpond/tximeta_GRCh38_ensembl/{coldata_id}/transcripts.rds",
        genes              = "out/fishpond/tximeta_GRCh38_ensembl/{coldata_id}/genes.rds",
        transcripts_counts = "out/fishpond/tximeta_GRCh38_ensembl/{coldata_id}/transcripts_counts.csv.gz",
        genes_counts       = "out/fishpond/tximeta_GRCh38_ensembl/{coldata_id}/genes_counts.csv.gz"
    params:
        outdir   = "out/fishpond/tximeta_GRCh38_ensembl/{coldata_id}",
        bfc_path = "out/TximetaBFC"
    wildcard_constraints:
        coldata_id = r"[-a-zA-Z0-9_]+"
    conda:
        "../envs/r_fishpond.yaml"
    shell:
        "Rscript {input.script} {input.coldata} {params.outdir} {params.bfc_path}"
