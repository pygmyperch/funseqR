# EXPORTED

#' Database connection utilities for funseqR
#'
#' These functions manage connections to the SQLite database that stores funseqR data.
#'

#' Create a new funseqR database
#'
#' This function creates a new SQLite database for storing funseqR data.
#' It sets up the required schema and initializes the database.
#'
#' @param db_path Character string specifying the path to the database file.
#'   If the file already exists, it will not be overwritten unless force=TRUE.
#' @param force Logical. If TRUE, the existing database will be overwritten. Default is FALSE.
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE.
#'
#' @return A connection object to the newly created database.
#'
#' @importFrom RSQLite SQLite
#' @importFrom DBI dbConnect dbDisconnect dbExecute
#'
#' @export
create_funseq_db <- function(db_path, force = FALSE, verbose = TRUE) {
  # Check if file exists
  if (file.exists(db_path) && !force) {
    stop("Database file already exists. Use force=TRUE to overwrite.")
  }

  # Delete existing file if force=TRUE
  if (file.exists(db_path) && force) {
    file.remove(db_path)
    if (verbose) message("Existing database file removed.")
  }

  # Create directory if it doesn't exist
  db_dir <- dirname(db_path)
  if (!dir.exists(db_dir)) {
    dir.create(db_dir, recursive = TRUE)
    if (verbose) message("Created directory: ", db_dir)
  }

  # Connect to database
  if (verbose) message("Creating new database: ", db_path)
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)

  # Initialize schema
  if (verbose) message("Initializing database schema...")
  create_funseq_schema(con, verbose = verbose)

  # Add package version info
  package_version <- as.character(packageVersion("funseqR"))
  r_version <- paste(R.version$major, R.version$minor, sep = ".")

  query <- paste0(
    "INSERT INTO metadata (key, value) VALUES ",
    "('funseqR_version', '", package_version, "'), ",
    "('r_version', '", r_version, "'), ",
    "('creation_date', '", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "')"
  )

  DBI::dbExecute(con, query)

  if (verbose) message("Database initialization complete.")

  return(con)
}

#' Connect to an existing funseqR database
#'
#' This function connects to an existing funseqR database.
#'
#' @param db_path Character string specifying the path to the database file.
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE.
#'
#' @return A connection object to the database.
#'
#' @importFrom RSQLite SQLite
#' @importFrom DBI dbConnect dbGetQuery
#'
#' @export
connect_funseq_db <- function(db_path, verbose = TRUE) {
  # Check if file exists
  if (!file.exists(db_path)) {
    stop("Database file does not exist: ", db_path)
  }

  # Connect to database
  if (verbose) message("Connecting to database: ", db_path)
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)

  # Check if it's a valid funseqR database
  tables <- DBI::dbListTables(con)
  required_tables <- c("metadata", "input_files")

  if (!all(required_tables %in% tables)) {
    DBI::dbDisconnect(con)
    stop("The file does not appear to be a valid funseqR database.")
  }

  # Get and display database info
  if (verbose) {
    metadata <- DBI::dbGetQuery(con, "SELECT key, value FROM metadata")
    version <- metadata$value[metadata$key == "funseqR_version"]
    creation_date <- metadata$value[metadata$key == "creation_date"]

    message("Connected to funseqR database (version: ", version, ")")
    message("Creation date: ", creation_date)
  }

  return(con)
}

#' Check if a database connection is open
#'
#' @param con A database connection object.
#'
#' @return Logical. TRUE if the connection is open, FALSE otherwise.
#'
#' @importFrom DBI dbIsValid
#'
#' @export
is_db_connected <- function(con) {
  return(DBI::dbIsValid(con))
}

#' Close a database connection
#'
#' @param con A database connection object.
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE.
#'
#' @return Logical. TRUE if the connection was successfully closed, FALSE otherwise.
#'
#' @importFrom DBI dbDisconnect
#'
#' @export
close_funseq_db <- function(con, verbose = TRUE) {
  if (DBI::dbIsValid(con)) {
    DBI::dbDisconnect(con)
    if (verbose) message("Database connection closed.")
    return(TRUE)
  } else {
    if (verbose) message("Connection already closed.")
    return(FALSE)
  }
}

#' Get the path of a database from its connection
#'
#' @param con A database connection object.
#'
#' @return Character string. The path to the database file.
#'
#' @importFrom DBI dbGetInfo
#'
#' @export
get_db_path <- function(con) {
  if (!DBI::dbIsValid(con)) {
    stop("Invalid database connection.")
  }

  return(DBI::dbGetInfo(con)$dbname)
}


# INTERNAL

#' Begin a transaction in the database
#'
#' @param con A database connection object.
#'
#' @return Invisible NULL.
#'
#' @importFrom DBI dbExecute
#'
begin_transaction <- function(con) {
  DBI::dbExecute(con, "BEGIN TRANSACTION")
  return(invisible(NULL))
}

#' Commit a transaction in the database
#'
#' @param con A database connection object.
#'
#' @return Invisible NULL.
#'
#' @importFrom DBI dbExecute
#'
commit_transaction <- function(con) {
  DBI::dbExecute(con, "COMMIT")
  return(invisible(NULL))
}

#' Rollback a transaction in the database
#'
#' @param con A database connection object.
#'
#' @return Invisible NULL.
#'
#' @importFrom DBI dbExecute
#'
rollback_transaction <- function(con) {
  DBI::dbExecute(con, "ROLLBACK")
  return(invisible(NULL))
}
