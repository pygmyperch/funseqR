#' Test connection to UniProt REST API
#'
#' Performs connectivity tests for the UniProt REST API by attempting a basic
#' connection and a simple query. Prints the results of response times, status
#' codes, and connection information to the console.
#'
#' @return This function does not return any object. It prints the test results
#'   to the console, including:
#'   \itemize{
#'     \item Basic connectivity test results
#'     \item Sample query test results
#'     \item Connection information (proxy settings, internet status)
#'   }
#' @importFrom httr GET status_code
#' @importFrom curl ie_get_proxy_for_url has_internet
#' @examples
#' \dontrun{
#' # Test UniProt connection and print results
#' test_uniprot_connection()
#' }
test_uniprot_connection <- function() {
  results <- list()
  
  # test basic connectivity
  message("Testing basic connectivity to UniProt...")
  
  tryCatch({
    start_time <- Sys.time()
    response <- httr::GET("https://rest.uniprot.org")
    end_time <- Sys.time()
    
    results$basic_connection <- list(
      success = TRUE,
      status = httr::status_code(response),
      response_time = as.numeric(difftime(end_time, start_time, units = "secs")),
      message = "Successfully connected to UniProt"
    )
  }, error = function(e) {
    results$basic_connection <- list(
      success = FALSE,
      status = NA,
      response_time = NA,
      message = paste("Connection error:", e$message)
    )
  })
  
  # test simple query
  message("Testing simple query...")
  
  tryCatch({
    start_time <- Sys.time()
    response <- httr::GET(
      url = "https://rest.uniprot.org/uniprotkb/search",
      query = list(
        query = "accession:P02879",
        fields = "accession,id",
        format = "tsv"
      )
    )
    end_time <- Sys.time()
    
    results$query_test <- list(
      success = httr::status_code(response) == 200,
      status = httr::status_code(response),
      response_time = as.numeric(difftime(end_time, start_time, units = "secs")),
      message = if(httr::status_code(response) == 200) "Query successful" else "Query failed"
    )
  }, error = function(e) {
    results$query_test <- list(
      success = FALSE,
      status = NA,
      response_time = NA,
      message = paste("Query error:", e$message)
    )
  })
  
  # print results
  message("\nConnection test results:")
  message(sprintf("Basic connection: %s (%.2f seconds)", 
                  ifelse(results$basic_connection$success, "SUCCESS", "FAILED"),
                  results$basic_connection$response_time))
  message(sprintf("Status code: %s", results$basic_connection$status))
  message(sprintf("Message: %s", results$basic_connection$message))
  
  message("\nQuery test results:")
  message(sprintf("Query test: %s (%.2f seconds)", 
                  ifelse(results$query_test$success, "SUCCESS", "FAILED"),
                  results$query_test$response_time))
  message(sprintf("Status code: %s", results$query_test$status))
  message(sprintf("Message: %s", results$query_test$message))
  
  # show connection info
  message("\nConnection information:")
  message("Proxy settings:")
  print(curl::ie_get_proxy_for_url("https://rest.uniprot.org"))
  message("Internet connection:", curl::has_internet())
  
  #return(results)
}

