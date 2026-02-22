# fishpond_tximeta.R
# Usage: Rscript fishpond_tximeta.R <coldata_tsv> <outdir> [bfc_path]
#
# coldata_tsv : TSV file with columns 'name' (sample label) and 'file'
#               (path to quant.sf produced by salmon with Gibbs samples)
# outdir      : output directory; receives transcripts.rds and genes.rds
# bfc_path    : (optional) BiocFileCache directory for tximeta reference
#               downloads; default: out/TximetaBFC

suppressPackageStartupMessages({
    library(tximeta)
    library(fishpond)
    library(org.Hs.eg.db)
    library(SummarizedExperiment)
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript fishpond_tximeta.R <coldata_tsv> <outdir> [bfc_path]")
}

coldata_tsv <- args[1]
outdir      <- args[2]
bfc_path    <- if (length(args) >= 3) args[3] else "out/TximetaBFC"

dir.create(outdir,   showWarnings = FALSE, recursive = TRUE)
dir.create(bfc_path, showWarnings = FALSE, recursive = TRUE)

setTximetaBFC(bfc_path)

coldata <- as.data.frame(fread(coldata_tsv, sep = "\t"))

if (!all(c("name", "file") %in% names(coldata))) {
    stop("coldata TSV must contain 'name' and 'file' columns")
}

# tximeta expects 'names' and 'files' columns
coldata$names <- coldata$name
coldata$files <- coldata$file

missing <- !file.exists(coldata$files)
if (any(missing)) {
    warning(
        sum(missing), " quant.sf file(s) not found and will be skipped:\n",
        paste(coldata$files[missing], collapse = "\n")
    )
    coldata <- coldata[!missing, ]
}

# Transcript-level SummarizedExperiment (GRCh38 Ensembl; release detected
# automatically from the salmon index metadata embedded in quant.sf)
se <- tximeta(coldata)

# Add gene-level SYMBOL annotations
se <- addIds(se, "SYMBOL", gene = TRUE)

# Scale inferential replicates (Gibbs samples) prior to downstream testing
se <- scaleInfReps(se)

saveRDS(se, file.path(outdir, "transcripts.rds"))
message("Saved: ", file.path(outdir, "transcripts.rds"))

# Gene-level summary
gse <- summarizeToGene(se)

saveRDS(gse, file.path(outdir, "genes.rds"))
message("Saved: ", file.path(outdir, "genes.rds"))
