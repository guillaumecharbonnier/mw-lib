# These functions define default arguments for Rmd templates defined in this folder.

knitPostDeseq2Template <- function(
                               d = vsd,
                               r = res,
                               p_cutoff = 0.05,
                               FC_cutoff = 2,
                               output_dir = book_from_rmd,
                               chunk_label_prefix = opts_current$get("label")) {
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
                               chunk_label_prefix = opts_current$get("label")) {
  src <- knitr::knit_expand("templates/wgcna.Rmd")
  res <- knitr::knit_child(
    text = unlist(src),
    envir = environment(),
    quiet = TRUE
  )
  cat(res, sep = "\n")
}