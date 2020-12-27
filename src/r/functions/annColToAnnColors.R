# This function solves limitations in default aheatmap() palette for annCol.
# Provide as input the data.frame you provide for aheatmap annCol arg
# and it return the list object you should provide to aheatmap annColors.
annColToAnnColors <- function(ann_col) {
  ann_colors <- list()
  seed_pal <- brewer.pal(n = 9, name = "Set1")

  pal_start <- 0
  pal_continuous <- 1
  for (col_name in colnames(ann_col)) {
    if (is.numeric(ann_col[[col_name]])){
        ann_colors[[col_name]] <- c(
          "#FFFFFF",
          seed_pal[[pal_continuous]]
        )
        pal_continuous <- (pal_continuous %% 9) + 1
      } else {
      pal <- createPalette(
        length(unique(ann_col[, col_name])),
        seed_pal[pal_start + 1]
      )
      names(pal) <- NULL
      ann_colors[[col_name]] <- pal
      pal_start <- (pal_start + 1) %% 9
    }
  }
  return(ann_colors)
}
