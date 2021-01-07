loadLibrary <- function(package) {
  if (!require(basename(package), character.only = TRUE)) BiocManager::install(package, update = FALSE)
  library(basename(package), character.only = TRUE)
}

loadLibraries <- function(packages) {
  invisible(lapply(packages, loadLibrary))
}
