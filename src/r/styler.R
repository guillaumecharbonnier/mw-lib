#!/usr/bin/env Rscript

styler::style_dir("functions",
  recursive = FALSE,
  filetype = c("R", "Rmd")
)
