
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

exportIntegrationDeaToXlsx <- function(d, filepath) {
  dir.create(
    dirname(filepath),
    recursive = TRUE,
    showWarnings = FALSE
  )
  # Export to xlsx
  absMax <- function(x){return(x[which.max(abs(x))])}
  d <- d[!is.na(d$hugo_symbol),]
  d_log2FoldChange <- reshape2::dcast(
    data = d,
    formula = hugo_symbol ~ design,
    value.var = "log2FoldChange",
    fun.aggregate = absMax,
    fill = 0
  )

  d_pvalue <- reshape2::dcast(
    data = d,
    formula = hugo_symbol ~ design,
    value.var = "pvalue",
    fun.aggregate = min
  )

  wb <- createWorkbook()
  addWorksheet(
    wb,
    "log2FoldChange"
  )
  addWorksheet(
    wb,
    "pvalue"
  )
  writeDataTable(
    wb,
    "log2FoldChange",
    d_log2FoldChange,
    startCol = 1
  )
  writeDataTable(
    wb,
    "pvalue",
    d_pvalue,
    startCol = 1
  )
  conditionalFormatting(
    wb,
    sheet = "pvalue",
    rows = 1:nrow(d_pvalue) + 1,
    type = "colourScale",
    style = c("darkgreen", "white"),
    rule = c(0, 0.1),
    cols = 2:ncol(d_pvalue)
  )
  conditionalFormatting(
    wb,
    sheet = "log2FoldChange",
    rows = 1:nrow(d_log2FoldChange) + 1,
    type = "colourScale",
    style = c("blue", "white", "red"),
    rule = c(-3, 0, 3),
    cols = 2:ncol(d_log2FoldChange)
  )
  for (worksheet in c("pvalue", "log2FoldChange")){
    freezePane(
      wb,
      sheet = worksheet,
      firstRow = TRUE,
      firstCol = TRUE
    )
    headerStyle <- createStyle(
      # fontSize = 18,
      # fontName = "Arial",
      # textDecoration = "bold",
      # halign = "left",
      # fgFill = "#1A33CC",
      # border = "TopBottomLeftRight",
      textRotation = 90
    )
    addStyle(
      wb,
      sheet = worksheet,
      style = headerStyle,
      rows = 1,
      cols = 1:ncol(d_pvalue)
    )
    setRowHeights(
      wb,
      sheet = worksheet,
      rows = 1,
      heights = 200
    )
    setColWidths(
      wb,
      sheet = worksheet,
      cols = 2:ncol(d_pvalue),
      widths = 4,
    )
    setColWidths(
      wb,
      sheet = worksheet,
      cols = 1,
      widths = 8,
    )
    saveWorkbook(
      wb,
      file = filepath,
      overwrite = TRUE
    )
  }
}