#' Extract BLAST database metadata using blastdbcmd
#'
#' This function uses the blastdbcmd utility to extract metadata about a BLAST database
#'
#' @param db_path Character string specifying the path to the BLAST database directory
#' @param db_name Character string specifying the name of the BLAST database
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE
#'
#' @return A list containing database metadata, or NULL if extraction fails
#'
#' @export
extract_blast_db_metadata <- function(db_path, db_name, verbose = TRUE) {
  # Construct full database path
  db_file <- file.path(db_path, db_name)

  if (verbose) message("Extracting metadata for database: ", db_file)

  # Check if blastdbcmd is available
  blastdbcmd_check <- system("which blastdbcmd", ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (blastdbcmd_check != 0) {
    if (verbose) message("blastdbcmd not found in PATH. Cannot extract database metadata.")
    return(NULL)
  }

  # Run blastdbcmd -info to get database information
  cmd <- paste("blastdbcmd -info -db", shQuote(db_file))
  if (verbose) message("Running command: ", cmd)

  tryCatch({
    # Capture output
    output <- system(cmd, intern = TRUE, ignore.stderr = FALSE)

    if (length(output) == 0) {
      if (verbose) message("No output from blastdbcmd -info")
      return(NULL)
    }

    # Parse the output
    metadata <- list(
      db_path = db_path,
      db_name = db_name,
      db_full_path = db_file,
      extraction_date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      raw_output = paste(output, collapse = "\n")
    )

    # Extract specific information from the output
    for (line in output) {
      # Database type
      if (grepl("Database:", line)) {
        metadata$db_title <- gsub("^Database:\\s*", "", line)
      }

      # Number of sequences
      if (grepl("sequences;", line)) {
        seq_match <- regexpr("[0-9,]+\\s+sequences", line)
        if (seq_match > 0) {
          seq_text <- regmatches(line, seq_match)
          metadata$num_sequences <- as.numeric(gsub("[^0-9]", "", seq_text))
        }
      }

      # Total length
      if (grepl("total length", line)) {
        len_match <- regexpr("[0-9,]+\\s+total length", line)
        if (len_match > 0) {
          len_text <- regmatches(line, len_match)
          metadata$total_length <- as.numeric(gsub("[^0-9]", "", len_text))
        }
      }

      # Date created (FIXED - handle case where date and longest sequence are on same line)
      if (grepl("Date:", line)) {
        # Extract just the date part, stopping at tab or "Longest sequence"
        date_part <- gsub("^Date:\\s*", "", line)
        # Split on tab or "Longest sequence" and take the first part
        date_part <- strsplit(date_part, "\\t|Longest sequence")[[1]][1]
        # Trim whitespace
        metadata$db_date <- trimws(date_part)
      }

      # Longest sequence (IMPROVED - handle case where it's on same line as date)
      if (grepl("Longest sequence:", line)) {
        # Extract the number before "residues"
        longest_match <- regexpr("[0-9,]+\\s+residues", line)
        if (longest_match > 0) {
          longest_text <- regmatches(line, longest_match)
          metadata$longest_sequence <- as.numeric(gsub("[^0-9]", "", longest_text))
        }
      }
    }

    # Try to extract version information from database files
    tryCatch({
      # Look for .nal or .pal files that might contain version info
      db_dir <- dirname(db_file)
      version_files <- list.files(db_dir, pattern = paste0(basename(db_name), "\\.(nal|pal)$"), full.names = TRUE)

      if (length(version_files) > 0) {
        # Try to read the first version file
        version_content <- readLines(version_files[1], n = 10, warn = FALSE)
        version_info <- paste(version_content, collapse = " ")

        # Look for version patterns
        version_match <- regexpr("version\\s+[0-9\\.]+", version_info, ignore.case = TRUE)
        if (version_match > 0) {
          metadata$db_version <- regmatches(version_info, version_match)
        }
      }
    }, error = function(e) {
      if (verbose) message("Could not extract version information: ", e$message)
    })

    # Ensure all fields are single values (not vectors)
    metadata <- lapply(metadata, function(x) {
      if (length(x) > 1) {
        paste(x, collapse = " ")
      } else if (length(x) == 0 || is.null(x)) {
        NA_character_
      } else {
        as.character(x)
      }
    })

    if (verbose) {
      message("Extracted database metadata:")
      message("  Title: ", metadata$db_title %||% "Unknown")
      message("  Sequences: ", format(as.numeric(metadata$num_sequences %||% 0), big.mark = ","))
      message("  Total length: ", format(as.numeric(metadata$total_length %||% 0), big.mark = ","))
      message("  Date: ", metadata$db_date %||% "Unknown")
      message("  Longest sequence: ", format(as.numeric(metadata$longest_sequence %||% 0), big.mark = ","))
    }

    return(metadata)

  }, error = function(e) {
    if (verbose) message("Error extracting database metadata: ", e$message)
    return(NULL)
  })
}

#' Store BLAST database metadata in the database
#'
#' @param con A database connection object
#' @param blast_param_id The ID of the BLAST parameters to associate with this metadata
#' @param metadata A list containing database metadata as returned by extract_blast_db_metadata
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE
#'
#' @return The ID of the stored metadata record, or NULL if storage fails
#'
#' @export
store_blast_db_metadata <- function(con, blast_param_id, metadata, verbose = TRUE) {
  if (is.null(metadata) || !is.list(metadata)) {
    if (verbose) message("No metadata to store")
    return(NULL)
  }

  # Ensure the metadata table exists
  ensure_blast_db_metadata_table(con, verbose = FALSE)

  # Check if metadata already exists for this blast_param_id
  existing <- DBI::dbGetQuery(
    con,
    "SELECT metadata_id FROM blast_database_metadata WHERE blast_param_id = ?",
    params = list(blast_param_id)
  )

  if (nrow(existing) > 0) {
    if (verbose) message("Database metadata already exists for blast_param_id ", blast_param_id)
    return(existing$metadata_id[1])
  }

  # Insert metadata
  tryCatch({
    DBI::dbExecute(
      con,
      "INSERT INTO blast_database_metadata (
        blast_param_id, db_path, db_name, db_full_path, db_title,
        num_sequences, total_length, db_date, db_version, longest_sequence,
        extraction_date, raw_output
      ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
      params = list(
        blast_param_id,
        metadata$db_path,
        metadata$db_name,
        metadata$db_full_path,
        metadata$db_title,
        metadata$num_sequences,
        metadata$total_length,
        metadata$db_date,
        metadata$db_version,
        metadata$longest_sequence,
        metadata$extraction_date,
        metadata$raw_output
      )
    )

    # Get the inserted record ID
    metadata_id <- DBI::dbGetQuery(
      con,
      "SELECT metadata_id FROM blast_database_metadata WHERE blast_param_id = ?",
      params = list(blast_param_id)
    )$metadata_id[1]

    if (verbose) message("Stored database metadata with ID: ", metadata_id)
    return(metadata_id)

  }, error = function(e) {
    if (verbose) message("Error storing database metadata: ", e$message)
    return(NULL)
  })
}

#' Ensure the blast_database_metadata table exists
#'
#' @param con A database connection object
#' @param verbose Logical. If TRUE, print progress information. Default is FALSE
#'
#' @return Invisible NULL
#'
#' @export
ensure_blast_db_metadata_table <- function(con, verbose = FALSE) {
  # Check if table exists
  tables <- DBI::dbListTables(con)
  if ("blast_database_metadata" %in% tables) {
    if (verbose) message("blast_database_metadata table already exists")
    return(invisible(NULL))
  }

  # Create table
  if (verbose) message("Creating blast_database_metadata table...")

  DBI::dbExecute(con, "
    CREATE TABLE blast_database_metadata (
      metadata_id INTEGER PRIMARY KEY,
      blast_param_id INTEGER NOT NULL,
      db_path TEXT NOT NULL,
      db_name TEXT NOT NULL,
      db_full_path TEXT NOT NULL,
      db_title TEXT,
      num_sequences INTEGER,
      total_length INTEGER,
      db_date TEXT,
      db_version TEXT,
      longest_sequence INTEGER,
      extraction_date TEXT NOT NULL,
      raw_output TEXT,
      FOREIGN KEY (blast_param_id) REFERENCES blast_parameters (blast_param_id)
    )
  ")

  DBI::dbExecute(con, "CREATE INDEX idx_blast_db_meta_param ON blast_database_metadata (blast_param_id)")
  DBI::dbExecute(con, "CREATE INDEX idx_blast_db_meta_name ON blast_database_metadata (db_name)")

  if (verbose) message("blast_database_metadata table created successfully")
  return(invisible(NULL))
}

#' Get BLAST database metadata from the database
#'
#' @param con A database connection object
#' @param blast_param_id Optional. The ID of the BLAST parameters. If NULL, returns all metadata
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE
#'
#' @return A data frame containing database metadata
#'
#' @export
get_blast_db_metadata <- function(con, blast_param_id = NULL, verbose = TRUE) {
  # Check if table exists
  tables <- DBI::dbListTables(con)
  if (!"blast_database_metadata" %in% tables) {
    if (verbose) message("blast_database_metadata table does not exist")
    return(data.frame())
  }

  # Build query
  if (is.null(blast_param_id)) {
    query <- "SELECT * FROM blast_database_metadata ORDER BY extraction_date DESC"
    params <- list()
  } else {
    query <- "SELECT * FROM blast_database_metadata WHERE blast_param_id = ?"
    params <- list(blast_param_id)
  }

  # Execute query
  tryCatch({
    result <- DBI::dbGetQuery(con, query, params = params)
    if (verbose && nrow(result) > 0) {
      message("Retrieved ", nrow(result), " database metadata record(s)")
    }
    return(result)
  }, error = function(e) {
    if (verbose) message("Error retrieving database metadata: ", e$message)
    return(data.frame())
  })
}

