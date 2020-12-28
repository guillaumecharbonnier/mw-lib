## A wrapper to saveWidget which compensates for arguable BUG in
## saveWidget which requires `file` to be in current working
## directory.
saveWidgetFix <- function(widget, file, ...) {

  wd <- getwd()
  on.exit(setwd(wd))
  outDir <- dirname(file)
  file <- basename(file)
  setwd(outDir)
  saveWidget(widget, file = file, ...)
}