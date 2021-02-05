library(data.table)
library(edgeR)

d <- fread(snakemake@input[["tsv"]])

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

m <- rpkm(d[,-1])
d <- data.table(
  d[,1],
  m
)

fwrite(
  x = d,
  file = snakemake@output[["rpkm"]],
  sep = "\t"
)
