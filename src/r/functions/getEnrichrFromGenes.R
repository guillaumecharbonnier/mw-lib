# getEnrichrUrlFromGenes <- function(genes) {
#   if (length(genes) == 0) {
#     return("No gene in this query to Enrichr")
#   } else {
#     req.body <- list(
#       list = paste(
#         genes,
#         collapse = "\n"
#       )
#     )
#     post.req <- tryCatch({
#       httr::POST(
#         "http://maayanlab.cloud/Enrichr/addList",
#         encode = "multipart",
#         body = I(req.body)
#       )
#     }, error = function(e) {
#       success <- FALSE
#       attempts <- 0
#       while (!success && attempts < 3) {
#         Sys.sleep(2)
#         post.req <- tryCatch({
#           httr::POST(
#             "http://maayanlab.cloud/Enrichr/addList",
#             encode = "multipart",
#             body = I(req.body)
#           )
#         }, error = function(e) {
#           success <<- FALSE
#         })
#         success <- !is.null(post.req)
#         attempts <- attempts + 1
#       }
#       post.req
#     })
#     url <- paste0(
#       "https://maayanlab.cloud/Enrichr/enrich?dataset=",
#       jsonlite::fromJSON(httr::content(
#         post.req,
#         "text"
#       ))$shortId
#     )
#   }
# }

# getEnrichrUrlFromGenes <- function(genes) {
#   if (length(genes) == 0) {
#     return("No gene in this query to Enrichr")
#   } else {
#     req.body <- list(
#       list = paste(
#         genes,
#         collapse = "\n"
#       )
#     )
#     attempts <- 1
#     success <- FALSE
#     while (!success && attempts < 3) {
#       post.req <- tryCatch(
#         {
#           httr::POST(
#             "http://maayanlab.cloud/Enrichr/addList",
#             encode = "multipart",
#             body = I(req.body)
#           )
#         },
#         error = function(e) {
#           success <<- FALSE
#         }
#       )
#       success <- !is.null(post.req)
#       if (httr::status_code(post.req) == 429) {success <- FALSE}
#       Sys.sleep(attempts * 2)
#       attempts <- attempts + 1
#     }
#     url <- paste0(
#       "https://maayanlab.cloud/Enrichr/enrich?dataset=",
#       jsonlite::fromJSON(httr::content(
#         post.req,
#         "text"
#       ))$shortId
#     )
#     return(url)
#   }
# }




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
    success <- FALSE
    attempts <- 1
    while (!success && attempts <= 3) {
      post.req <- tryCatch(
        {
          httr::POST(
            "http://maayanlab.cloud/Enrichr/addList",
            encode = "multipart",
            body = I(req.body)
          )
        },
        error = function(e) {
          success <<- FALSE
        }
      )
      success <- !is.null(post.req)
      if (httr::status_code(post.req) == 429) {success <- FALSE}
      Sys.sleep(attempts * 2)
      attempts <- attempts + 1
    }
    if (!success) {
      return("Query to EnrichR failed 3 times, please try again later")
    } else {
      url <- paste0(
        "https://maayanlab.cloud/Enrichr/enrich?dataset=",
        jsonlite::fromJSON(httr::content(
          post.req,
          "text"
        ))$shortId
      )
      return(url)
    }
  }
}