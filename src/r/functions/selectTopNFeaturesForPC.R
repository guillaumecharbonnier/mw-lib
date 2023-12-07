# A function to select the top N features for a given PC
# sign can be either "positive" or "negative"
# d is the PCA object
# subset is an index vector to subset the features of d. If NULL, all features are used.
# this was designed to retrieve the top N features for a given dataset, eg with subset = rowData(d_pca)$dataset %in% "RNA"
selectTopNFeaturesForPC <- function(
  d = d_pca,
  subset = NULL,
  pc = 1,
  sign = c("positive", "negative"),
  n = 200
) {
  if (!is.null(subset)) {
    rot <- d$rotation[subset,pc]
  } else {
    rot <- d$rotation[,pc]
  }
  if (sign == "positive") {
    ranks <- rank(-rot)
  } else if (sign == "negative") {
    ranks <- rank(rot)
  } else {
    stop("sign must be either 'positive' or 'negative'")
  }
  ranks <- ranks[order(ranks)]

  return(names(ranks)[1:n])
}
# selectTopNFeaturesForPC(sign = "positive")
# selectTopNFeaturesForPC(
#   d = d_pca,
#   pc = 2,
#   sign = "negative",
#   n = 100
# )