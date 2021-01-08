#!/usr/bin/env Rscript

styler::style_dir(
  path = "functions",
  recursive = FALSE,
  filetype = c("R", "Rmd")
)

styler::style_dir(
  path = "../rmd/templates",
  recursive = FALSE,
  filetype = c("R", "Rmd")
)
