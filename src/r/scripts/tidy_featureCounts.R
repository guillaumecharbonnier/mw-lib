library(data.table)
library(edgeR)

d <- fread(snakemake@input[["tsv"]])

gene_length <- d$Length

d[, c("Chr", "Start", "End", "Strand", "Length") := NULL]

setnames(
  d,
  old = names(d[,-1]),
  new = gsub(
    pattern = "^.*/|.bam$",
    replacement = "",
    x = names(d[,-1])
  )
)

fwrite(
  x = d,
  file = snakemake@output[["raw"]],
  sep = "\t"
)

m <- rpkm(d[,-1], gene_length)
d_rpkm <- data.table(
  d[,1],
  m
)

fwrite(
  x = d_rpkm,
  file = snakemake@output[["rpkm"]],
  sep = "\t"
)


rpk <- as.matrix(d[,-1] / (gene_length / 1000))
dge_rpk <- DGEList(counts = rpk)
dge_rpk <- calcNormFactors(
  dge_rpk,
  method = "TMM"
)
tmm_rpk <- cpm(dge_rpk)
d_tmm_rpk <- data.table(
  Geneid = d$Geneid,
  tmm_rpk
)
fwrite(
  x = d_tmm_rpk,
  file = snakemake@output[["tmm_rpk"]],
  sep = "\t"
)
