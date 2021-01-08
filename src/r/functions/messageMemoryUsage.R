# Message memory usage (to dev and debug ram consmption)
messageMemoryUsage <- function() {
  x <- sort(
    sapply(
      ls(envir = parent.env(environment())),
      function(x) {
        object.size(get(x))
      }
    )
  )
  message(
    paste0(
      capture.output(x),
      collapse = "\n"
    )
  )
}
