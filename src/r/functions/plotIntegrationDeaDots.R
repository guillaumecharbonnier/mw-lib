
# pval_cutoff : + or - sign will be display for a gene-design pair only if the pvalue is below this cutoff
plotIntegrationDeaDots <- function(
  d,
  y_col = "hugo_symbol",
  x_col = "design",
  pvalue_cutoff = 0.05,
  log2FoldChange_cutoff = 0.1,
  squishOobValue = TRUE){
  d$signLog2FoldChange <- ifelse(
    d$log2FoldChange > log2FoldChange_cutoff &
      d$pvalue < pvalue_cutoff,
    "+",
    ifelse(
      d$log2FoldChange < - log2FoldChange_cutoff &
        d$pvalue < pvalue_cutoff,
      "–",
      ""
    )
  )
  d$mlog10pvalue <- ifelse(
    is.na(d$pvalue),
    0,
    ifelse(
      d$pvalue < 10^-10,
      10,
      -log10(d$pvalue)
    )
  )
  p <- ggplot(
    data.frame(d),
    aes(
      x = get(x_col),
      y = get(y_col),
      size = mlog10pvalue,
      label = signLog2FoldChange,
      color = log2FoldChange
    )
  )
  p <- p + geom_point()
  p <- p + geom_text(color="black")
  if (squishOobValue) {
    p <- p + scale_color_gradient2(
      low = "blue",
      mid = "lightgrey",
      high = "red",
      midpoint = 0,
      limits = c(-2,2),
      oob = scales::squish
    )
  } else {
    p <- p + scale_color_gradient2(
      low = "blue",
      mid = "lightgrey",
      high = "red",
      midpoint = 0
    )
  }
  p <- p + theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  )
  p <- p + xlab(x_col)
  p <- p + ylab(y_col)
  p
}
