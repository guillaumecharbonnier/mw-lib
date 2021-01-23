scaleRows <- function(x) {
  x <- t(scale(t(x)))
  # scaling can create NA rows when there is no variance
  x[is.na(x)] <- 0
  return(x)
}
