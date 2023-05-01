selectPcaKeyFeatures <- function(
  d = d_pca,
  n_pc_to_use = NULL,
  frac_of_var = 2) {
  if (is.null(n_pc_to_use)) {
    # Use the PC that account for approx half the variance in the dataset
    n_pc_to_use <- sum(cumsum(d$sdev) < sum(d$sdev) / frac_of_var)
  }
  rank_by_pc <- apply(
    d$rotation[,1:n_pc_to_use],
    2,
    rank
  )

  rank_by_pc <- lapply(
    seq_len(ncol(rank_by_pc)),
    function(i) rank_by_pc[,i]
  )

  res <- mapply(
    FUN = function(x,y) {x < y | x > max(x) - y},
    x = rank_by_pc,
    y = d$sdev[1:n_pc_to_use]
  )

  unique(rownames(res)[rowSums(res) > 0])
}