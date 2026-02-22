"""
Aim:
    Rules to import salmon quant.sf files with tximeta (human GRCh38 / Ensembl)
    and run fishpond preprocessing, producing:
      - transcript-level SummarizedExperiment (transcripts.rds)
      - gene-level SummarizedExperiment      (genes.rds)

Usage:
    Register a coldata TSV in mwconf['ids'] under a chosen <coldata_id> key.
    The TSV must have two tab-separated columns:
        name    sample label (used as column name in the SE)
        file    path to the corresponding salmon quant.sf

    Then request:
        out/fishpond/tximeta_GRCh38_ensembl/<coldata_id>/transcripts.rds
        out/fishpond/tximeta_GRCh38_ensembl/<coldata_id>/genes.rds
"""


def _fishpond_tximeta_quant_files(wildcards):
    """Return the list of quant.sf paths declared in the coldata TSV."""
    coldata_path = mwconf['ids'][wildcards.coldata_id]
    with open(coldata_path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        return [row['file'] for row in reader]


rule fishpond_tximeta_GRCh38_ensembl:
    """
    Created:
        2026-02-22
    Aim:
        Import salmon quant.sf files (human GRCh38, Ensembl release detected
        automatically from salmon index metadata) with tximeta, add SYMBOL IDs,
        scale inferential replicates, and summarise to gene level.
    Input coldata_id:
        Key in mwconf['ids'] that points to a TSV with 'name' and 'file' columns.
    Test:
        out/fishpond/tximeta_GRCh38_ensembl/my_coldata_id/transcripts.rds
        out/fishpond/tximeta_GRCh38_ensembl/my_coldata_id/genes.rds
    """
    input:
        coldata  = lambda wildcards: mwconf['ids'][wildcards.coldata_id],
        quant_sf = _fishpond_tximeta_quant_files,
        script   = "src/r/scripts/fishpond_tximeta.R"
    output:
        transcripts = "out/fishpond/tximeta_GRCh38_ensembl/{coldata_id}/transcripts.rds",
        genes       = "out/fishpond/tximeta_GRCh38_ensembl/{coldata_id}/genes.rds"
    params:
        outdir   = "out/fishpond/tximeta_GRCh38_ensembl/{coldata_id}",
        bfc_path = "out/TximetaBFC"
    wildcard_constraints:
        coldata_id = r"[-a-zA-Z0-9_]+"
    conda:
        "../envs/r_fishpond.yaml"
    shell:
        "Rscript {input.script} {input.coldata} {params.outdir} {params.bfc_path}"
