getEnrichrUrlFromGenes <- function(genes) {
  if (length(genes) == 0) {
    return("No gene in this query to Enrichr")
  } else {
    req.body <- list(
      list = paste(
        genes,
        collapse = "\n"
      )
    )
    post.req <- httr::POST(
      "http://maayanlab.cloud/Enrichr/addList",
      encode = "multipart",
      body = I(req.body)
    )
    url <- paste0(
      "https://maayanlab.cloud/Enrichr/enrich?dataset=",
      jsonlite::fromJSON(httr::content(
        post.req,
        "text"
      ))$shortId
    )
  }
}
