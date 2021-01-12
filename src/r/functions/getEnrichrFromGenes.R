getEnrichrUrlFromGenes <- function(genes) {
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
    "https://amp.pharm.mssm.edu/Enrichr/enrich?dataset=",
    jsonlite::fromJSON(httr::content(
      post.req,
      "text"
    ))$shortId
  )
}
