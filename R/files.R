# EXPORTED

#' File management functions for funseqR
#'
#' These functions manage input files within the funseqR database.
#'

#' Register an input file in the database
#'
#' This function adds information about an input file to the database.
#'
#' @param con A database connection object.
#' @param file_path Character string specifying the path to the file.
#' @param file_type Character string specifying the type of the file (e.g., "vcf", "fasta").
#' @param calculate_hash Logical. If TRUE, calculate and store a hash of the file. Default is TRUE.
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE.
#'
#' @return The ID of the newly registered file.
#'
#' @importFrom DBI dbExecute dbGetQuery
#' @export
register_input_file <- function(con, file_path, file_type, calculate_hash = TRUE, verbose = TRUE) {

  # Check if file exists
  if (!file.exists(file_path)) {
    stop("File does not exist: ", file_path)
  }

  # Get file name
  file_name <- basename(file_path)

  # Calculate hash if requested
  file_hash <- NULL
  if (calculate_hash) {
    if (verbose) message("Calculating file hash...")
    file_hash <- digest::digest(file_path, algo = "sha256", file = TRUE)
  }

  # Check if this file has already been registered
  existing <- DBI::dbGetQuery(
    con,
    "SELECT file_id FROM input_files
     WHERE file_name = ? AND file_path = ?",
    params = list(file_name, file_path)
  )

  if (nrow(existing) > 0) {
    if (verbose) message("File already registered with ID ", existing$file_id[1])
    return(existing$file_id[1])
  }

  # Register file
  current_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  DBI::dbExecute(
    con,
    "INSERT INTO input_files (file_type, file_name, file_path, file_hash, import_date)
     VALUES (?, ?, ?, ?, ?)",
    params = list(file_type, file_name, file_path, file_hash, current_time)
  )

  # Get the ID of the newly registered file
  file_id <- DBI::dbGetQuery(
    con,
    "SELECT file_id FROM input_files
     WHERE file_name = ? AND file_path = ?",
    params = list(file_name, file_path)
  )$file_id[1]

  if (verbose) message("Registered file '", file_name, "' with ID ", file_id)

  return(file_id)
}

#' List all input files in the database
#'
#' @param con A database connection object.
#'
#' @return A data frame containing information about all input files in the database.
#'
#' @importFrom DBI dbGetQuery
#' @export
list_input_files <- function(con) {
  # Get all files
  DBI::dbGetQuery(
    con,
    "SELECT * FROM input_files ORDER BY file_id"
  )
}

#' Get information about a specific file
#'
#' @param con A database connection object.
#' @param file_id The ID of the file.
#'
#' @return A data frame containing information about the file.
#'
#' @importFrom DBI dbGetQuery
#' @export
get_input_file <- function(con, file_id) {
  file_info <- DBI::dbGetQuery(
    con,
    "SELECT * FROM input_files WHERE file_id = ?",
    params = list(file_id)
  )

  if (nrow(file_info) == 0) {
    stop("File with ID ", file_id, " not found.")
  }

  return(file_info)
}