# Call from r_aggregate_tables.smk
# Take different tables and merge them based on their common columns
# Initially created to merge count tables from subread featureCounts

library(data.table)

files <- fread(snakemake@input)

# files <- c(
#   "out/subread/featureCounts_-O_-t_exon_-g_gene_name_gtf-GRCh38-ensembl_bam-GRCh38-exp-RNA-Thymoglobulin.tsv",
#   "out/subread/featureCounts_-O_-t_exon_-g_gene_name_gtf-GRCh38-ensembl_bam-GRCh38-exp-RNA-HHEX.tsv",
#   "out/subread/featureCounts_-O_-t_exon_-g_gene_name_gtf-GRCh38-ensembl_bam-GRCh38-exp-RNA-TALL-Agilent-XT-HS2.tsv",
#   "out/subread/featureCounts_-O_-t_exon_-g_gene_name_gtf-GRCh38-ensembl_bam-GRCh38-exp-RNA-Thymocyte.tsv"
# )

d <- list()
for (file in files) {
  d[[file]] <- fread(file)
}

dm <- Reduce(
  function(...) merge(
    ...,
    all=TRUE
  ),
  d
)
fwrite(
  x = dm,
  file = snakemake@output,
)
