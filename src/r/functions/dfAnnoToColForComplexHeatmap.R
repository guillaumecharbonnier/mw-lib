colorblind_palette_8 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

colorblind_palette_15 <- c(
  "#000000","#004949","#009292","#ff6db6","#ffb6db",
  "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
  "#920000","#924900","#db6d00","#24ff24","#ffff6d"
)

# This function is used to set colors for a discrete dataframe used as annotations for ComplexHeatmap
#
# Example:
# df_col <- dfAnnoToColForComplexHeatmap(df_anno)
# HeatmapAnnotation(
#   df = df_anno,
#   col = df_col
# )

vecAnnoToColForComplexHeatmap <- function(x) {
  if (is.numeric(x)) {
    return(NULL)
  } else {
    # For discrete variables, assign a unique color to each unique value
    y <- sort(unique(x))
    if (length(y) > 15) {
      warning("More than 15 different levels. Annotation heatmap will be unreadable")
    }
    if (length(y) > 8) {
      pal = colorblind_palette_15
    } else {
      pal = colorblind_palette_8
    }
    colors = structure(
      pal[1:length(y)],
      names = as.character(y)
    )
  }
  return(colors)
}

dfAnnoToColForComplexHeatmap <- function(
  df,
  pal = NULL
) {
  anno_col <- lapply(
    as.list(df),
    vecAnnoToColForComplexHeatmap
  )
  # we remove the NULL elements
  anno_col <- anno_col[sapply(anno_col, Negate(is.null))]
  return(anno_col)
}