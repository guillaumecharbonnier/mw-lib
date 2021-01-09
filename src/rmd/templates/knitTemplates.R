# These functions define default arguments for Rmd templates defined in this folder.


#' @params cols is a vector of column names to keep in the produced DataTable.
knitPostDeseq2Template <- function(
                                   d = vsd,
                                   r = res,
                                   p_cutoff = 0.05,
                                   FC_cutoff = 2,
                                   cols = c(
                                     "ensembl_gene_id",
                                     "hugo_symbol",
                                     "baseMean",
                                     "log2FoldChange",
                                     "pvalue",
                                     "padj"
                                   ),
                                   output_dir = book_from_rmd,
                                   chunk_label_prefix = opts_current$get("label")) {
  required_r_cols <- c(
    "baseMean",
    "log2FoldChange",
    "pvalue",
    "padj"
  )

  if (class(d) != "DESeqDataSet") {
    error("d should be a 'DESeqDataSet' object")
  }
  if (class(r) != "DESeqResults") {
    error("r should be a 'DESeqResults' object produced for example by lfcShrink()")
  }
  if (!"hugo_symbol" %in% names(rowData(d))) {
    warning("hugo_symbol column is missing. It can be useful to produce enrichment analyses.")
  }


  src <- knitr::knit_expand("templates/postdeseq2.Rmd")
  res <- knitr::knit_child(
    text = unlist(src),
    envir = environment(),
    quiet = TRUE
  )
  cat(res, sep = "\n")
}


knitWgcnaTemplate <- function(
                              d = vsd,
                              selected_power = NULL,
                              output_dir = book_from_rmd,
                              chunk_label_prefix = opts_current$get("label")) {
  src <- knitr::knit_expand("templates/wgcna.Rmd")
  res <- knitr::knit_child(
    text = unlist(src),
    envir = environment(),
    quiet = TRUE
  )
  cat(res, sep = "\n")
}
