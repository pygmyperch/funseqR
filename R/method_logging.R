# INTERNAL

#' Log method execution details for reproducibility
#'
#' This function logs the execution details of funseqR functions to enable
#' complete reproducibility and transparency. It captures command text,
#' parameters, execution time, and software environment information.
#'
#' @param con Database connection object
#' @param method_type Character. Type of method: "blast", "annotation", "enrichment", "sequence", "statistical", "import", "export"
#' @param function_name Character. Name of the function being executed
#' @param command_text Character. The actual command or operation being performed (optional)
#' @param parameters List. Parameters passed to the function (will be converted to JSON)
#' @param execution_time_seconds Numeric. Time taken to execute the function in seconds (optional)
#' @param success Logical. Whether the function executed successfully (default TRUE)
#' @param error_message Character. Error message if execution failed (optional)
#' @param verbose Logical. Print logging information (default FALSE)
#'
#' @return Invisible NULL
#'
#' @details
#' This function is called internally by funseqR functions to maintain a complete
#' audit trail of all operations performed. The logged information includes:
#' \itemize{
#'   \item Method type and function name
#'   \item Command text (e.g., BLAST command line, API URLs)
#'   \item Parameters as JSON
#'   \item Execution time and success status
#'   \item R version and package versions
#'   \item Timestamp of execution
#' }
#'
#' @importFrom DBI dbExecute
#' @importFrom jsonlite toJSON
#' @keywords internal
log_method_execution <- function(con, method_type, function_name, 
                                command_text = NULL, parameters = NULL,
                                execution_time_seconds = NULL, success = TRUE,
                                error_message = NULL, verbose = FALSE) {
  
  # Validate connection
  if (!DBI::dbIsValid(con)) {
    if (verbose) warning("Invalid database connection - method logging skipped")
    return(invisible(NULL))
  }
  
  # Check if method_log table exists
  existing_tables <- DBI::dbListTables(con)
  if (!"method_log" %in% existing_tables) {
    if (verbose) warning("method_log table not found - method logging skipped")
    return(invisible(NULL))
  }
  
  # Prepare parameters JSON
  parameters_json <- NULL
  if (!is.null(parameters)) {
    tryCatch({
      parameters_json <- jsonlite::toJSON(parameters, auto_unbox = TRUE, pretty = FALSE)
    }, error = function(e) {
      if (verbose) warning("Error converting parameters to JSON: ", e$message)
      parameters_json <- paste("Error converting parameters:", e$message)
    })
  }
  
  # Get R version and package versions
  r_version <- paste(R.version$major, R.version$minor, sep = ".")
  
  # Get key package versions
  package_versions <- tryCatch({
    key_packages <- c("funseqR", "DBI", "RSQLite", "httr", "jsonlite", "clusterProfiler", 
                      "vcfR", "Biostrings", "GenomicRanges")
    installed_packages <- utils::installed.packages()
    versions <- list()
    
    for (pkg in key_packages) {
      if (pkg %in% rownames(installed_packages)) {
        versions[[pkg]] <- installed_packages[pkg, "Version"]
      }
    }
    
    if (length(versions) > 0) {
      jsonlite::toJSON(versions, auto_unbox = TRUE, pretty = FALSE)
    } else {
      NULL
    }
  }, error = function(e) {
    if (verbose) warning("Error getting package versions: ", e$message)
    NULL
  })
  
  # Get current timestamp
  execution_date <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  # Insert log entry
  tryCatch({
    # Ensure all parameters are scalar values for DBI::dbExecute
    # Convert everything to simple types and handle NULL/vector issues
    safe_method_type <- as.character(method_type)[1]
    safe_function_name <- as.character(function_name)[1]
    safe_command_text <- if(is.null(command_text) || length(command_text) == 0) "" else as.character(command_text)[1]
    safe_parameters_json <- if(is.null(parameters_json) || length(parameters_json) == 0) "" else as.character(parameters_json)[1]
    safe_execution_date <- as.character(execution_date)[1]
    safe_r_version <- as.character(r_version)[1]
    safe_package_versions <- if(is.null(package_versions) || length(package_versions) == 0) "" else as.character(package_versions)[1]
    safe_execution_time <- if(is.null(execution_time_seconds) || length(execution_time_seconds) == 0) as.numeric(NA) else as.numeric(execution_time_seconds)[1]
    safe_success <- as.integer(as.logical(success)[1])  # Convert boolean to integer for SQLite
    safe_error_message <- if(is.null(error_message) || length(error_message) == 0) "" else as.character(error_message)[1]
    
    DBI::dbExecute(con, "
      INSERT INTO method_log (
        method_type, function_name, command_text, parameters_json,
        execution_date, r_version, package_versions, execution_time_seconds,
        success, error_message
      ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    ", list(
      safe_method_type, safe_function_name, safe_command_text, safe_parameters_json,
      safe_execution_date, safe_r_version, safe_package_versions, safe_execution_time,
      safe_success, safe_error_message
    ))
    
    if (verbose) {
      message("Method execution logged: ", method_type, " - ", function_name)
    }
  }, error = function(e) {
    if (verbose) warning("Error logging method execution: ", e$message)
  })
  
  return(invisible(NULL))
}

#' Wrapper function to time and log method execution
#'
#' This function wraps the execution of another function, automatically
#' timing it and logging the results.
#'
#' @param con Database connection object
#' @param method_type Character. Type of method being executed
#' @param function_name Character. Name of the function being executed
#' @param command_text Character. The actual command or operation (optional)
#' @param parameters List. Parameters passed to the function
#' @param func Function. The function to execute
#' @param verbose Logical. Print logging information (default FALSE)
#'
#' @return The result of the executed function
#'
#' @details
#' This wrapper function:
#' \itemize{
#'   \item Starts a timer
#'   \item Executes the provided function
#'   \item Captures execution time and success/failure
#'   \item Logs all information to the method_log table
#' }
#'
#' @keywords internal
execute_and_log <- function(con, method_type, function_name, 
                           command_text = NULL, parameters = NULL,
                           func, verbose = FALSE) {
  
  start_time <- Sys.time()
  success <- TRUE
  error_message <- NULL
  result <- NULL
  
  # Execute the function
  tryCatch({
    result <- func()
  }, error = function(e) {
    success <<- FALSE
    error_message <<- e$message
    if (verbose) warning("Function execution failed: ", e$message)
  })
  
  # Calculate execution time
  end_time <- Sys.time()
  execution_time_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Log the execution
  log_method_execution(
    con = con,
    method_type = method_type,
    function_name = function_name,
    command_text = command_text,
    parameters = parameters,
    execution_time_seconds = execution_time_seconds,
    success = success,
    error_message = error_message,
    verbose = verbose
  )
  
  # Re-throw error if function failed
  if (!success) {
    stop(error_message)
  }
  
  return(result)
}

#' Get summary of recent method executions
#'
#' Quick utility function to view recent method executions for debugging
#'
#' @param con Database connection object
#' @param limit Integer. Number of recent executions to retrieve (default 10)
#'
#' @return Data frame with recent method executions
#'
#' @keywords internal
get_recent_method_logs <- function(con, limit = 10) {
  
  if (!DBI::dbIsValid(con)) {
    stop("Invalid database connection")
  }
  
  # Check if method_log table exists
  existing_tables <- DBI::dbListTables(con)
  if (!"method_log" %in% existing_tables) {
    return(data.frame(message = "method_log table not found"))
  }
  
  # Get recent logs
  tryCatch({
    DBI::dbGetQuery(con, paste("
      SELECT 
        log_id, method_type, function_name, execution_date,
        execution_time_seconds, success, 
        substr(command_text, 1, 100) as command_preview
      FROM method_log
      ORDER BY execution_date DESC
      LIMIT", limit))
  }, error = function(e) {
    data.frame(message = paste("Error querying method logs:", e$message))
  })
}