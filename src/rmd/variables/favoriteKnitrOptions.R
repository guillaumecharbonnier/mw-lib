eval_after <- c(
  "fig.cap",
  "fig.height",
  "fig.width",
  "out.height",
  "out.width"
)

knitr::opts_chunk$set(
  echo = FALSE,
  message = FALSE,
  eval.after = eval_after,
  cache = TRUE
)
