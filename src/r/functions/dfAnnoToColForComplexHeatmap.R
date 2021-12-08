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
dfAnnoToColForComplexHeatmap <- function(df) {
  anno_col <- lapply(
    as.list(df),
    function(x) {
      y <- sort(unique(x))
      if (length(y) > 15) {
        warning("More than 15 different levels. Annotation heatmap will be unreadable")
      }
      colors = structure(
        colorblind_palette_15[1:length(y)],
        names = y
      )
    }
  )
  return(anno_col)
}