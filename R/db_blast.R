#' BLAST search and results management functions for funseqR
#'
#' These functions manage BLAST parameters and results within the funseqR database.
#'

#' Register BLAST parameters in the database
#'
#' This function stores BLAST search parameters in the database.
#'
#' @param con A database connection object.
#' @param project_id The ID of the project to associate with the BLAST search.
#' @param blast_type Character string specifying the type of BLAST search.
#'   Must be either "blastn" or "blastx".
#' @param db_name Character string specifying the name of the BLAST database.
#' @param db_path Character string specifying the path to the BLAST database.
#' @param e_value Numeric value specifying the E-value threshold for hits.
#' @param max_hits Integer specifying the maximum number of hits to return per query.
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE.
#'
#' @return The ID of the newly registered BLAST parameters.
#'
#' @importFrom DBI dbExecute dbGetQuery
#' @export
register_blast_params <- function(con, project_id, blast_type, db_name, db_path,
                           e_value, max_hits, verbose = TRUE) {
  # Check if project exists
  project <- DBI::dbGetQuery(
    con,
    "SELECT project_id FROM projects WHERE project_id = ?",
    params = list(project_id)
  )

  if (nrow(project) == 0) {
    stop("Project with ID ", project_id, " not found.")
  }

  # Validate blast_type
  blast_type <- match.arg(blast_type, c("blastn", "blastx"))

  # Register BLAST parameters
  current_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  DBI::dbExecute(
    con,
    "INSERT INTO blast_parameters (project_id, blast_type, db_name, db_path, e_value, max_hits, execution_date)
     VALUES (?, ?, ?, ?, ?, ?, ?)",
    params = list(project_id, blast_type, db_name, db_path, e_value, max_hits, current_time)
  )

  # Get the ID of the newly registered parameters
  param_id <- DBI::dbGetQuery(
    con,
    "SELECT blast_param_id FROM blast_parameters
     WHERE project_id = ? AND execution_date = ?
     ORDER BY blast_param_id DESC LIMIT 1",
    params = list(project_id, current_time)
  )$blast_param_id[1]

  if (verbose) message("Registered BLAST parameters with ID ", param_id)

  return(param_id)
}

#' Get BLAST parameters from the database
#'
#' @param con A database connection object.
#' @param blast_param_id The ID of the BLAST parameters to retrieve.
#'
#' @return A data frame containing the BLAST parameters.
#'
#' @importFrom DBI dbGetQuery
#' @export
get_blast_params <- function(con, blast_param_id) {
  params <- DBI::dbGetQuery(
    con,
    "SELECT * FROM blast_parameters WHERE blast_param_id = ?",
    params = list(blast_param_id)
  )

  if (nrow(params) == 0) {
    stop("BLAST parameters with ID ", blast_param_id, " not found.")
  }

  return(params)
}

#' List BLAST parameters for a project
#'
#' @param con A database connection object.
#' @param project_id The ID of the project.
#'
#' @return A data frame containing all BLAST parameters for the project.
#'
#' @importFrom DBI dbGetQuery
#' @export
list_blast_params <- function(con, project_id) {
  # Check if project exists
  project <- DBI::dbGetQuery(
    con,
    "SELECT project_id FROM projects WHERE project_id = ?",
    params = list(project_id)
  )

  if (nrow(project) == 0) {
    stop("Project with ID ", project_id, " not found.")
  }

  # Get parameters
  DBI::dbGetQuery(
    con,
    "SELECT * FROM blast_parameters WHERE project_id = ? ORDER BY blast_param_id",
    params = list(project_id)
  )
}

#' Import BLAST results into the database
#'
#' This function imports BLAST results from a tabular output file into the database.
#'
#' @param con A database connection object.
#' @param blast_param_id The ID of the BLAST parameters associated with these results.
#' @param results_file A character string specifying the path to the BLAST results file.
#' @param flanking_ids A vector of flanking sequence IDs that correspond to the query sequences.
#'   Must be in the same order as the queries in the BLAST results file.
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE.
#'
#' @return The number of BLAST results imported.
#'
#' @importFrom DBI dbExecute dbGetQuery
#' @importFrom progress progress_bar
#' @export
import_blast_results <- function(con, blast_param_id, results_file, flanking_ids, verbose = TRUE) {
  # Check if BLAST parameters exist
  params <- DBI::dbGetQuery(
    con,
    "SELECT * FROM blast_parameters WHERE blast_param_id = ?",
    params = list(blast_param_id)
  )

  if (nrow(params) == 0) {
    stop("BLAST parameters with ID ", blast_param_id, " not found.")
  }

  # Check if results file exists
  if (!file.exists(results_file)) {
    stop("BLAST results file does not exist: ", results_file)
  }

  # Read results file
  if (verbose) message("Reading BLAST results file...")
  blast_columns <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                    "qstart", "qend", "sstart", "send", "evalue", "bitscore")

  results <- read.table(results_file, col.names = blast_columns, stringsAsFactors = FALSE)

  if (nrow(results) == 0) {
    if (verbose) message("No BLAST hits found in the results file.")
    return(0)
  }

  # Match query IDs to flanking IDs
  if (verbose) message("Matching query IDs to flanking sequence IDs...")

  # Create a mapping from query IDs to flanking IDs
  # This assumes the query IDs in the BLAST results match the order of flanking_ids
  unique_queries <- unique(results$qseqid)

  if (length(unique_queries) > length(flanking_ids)) {
    stop("More query sequences in BLAST results than provided flanking IDs.")
  }

  # Map query IDs to flanking IDs
  query_to_flanking <- data.frame(
    qseqid = unique_queries,
    flanking_id = flanking_ids[1:length(unique_queries)],
    stringsAsFactors = FALSE
  )

  # Add flanking IDs to results
  results <- merge(results, query_to_flanking, by = "qseqid")

  # Start transaction
  DBI::dbExecute(con, "BEGIN TRANSACTION")

  # Set up progress bar
  if (verbose) {
    pb <- progress::progress_bar$new(
      format = "[:bar] :percent ETA: :eta",
      total = nrow(results),
      clear = FALSE,
      width = 60
    )
  }

  # Import results in batches
  batch_size <- 1000
  num_batches <- ceiling(nrow(results) / batch_size)

  tryCatch({
    for (i in 1:num_batches) {
      start_idx <- (i - 1) * batch_size + 1
      end_idx <- min(i * batch_size, nrow(results))
      batch <- results[start_idx:end_idx, ]

      # Extract description from sseqid if available
      hit_descriptions <- rep(NA_character_, nrow(batch))

      # For SwissProt format: "sp|P12345|GENE_SPECIES Description"
      sp_pattern <- "^sp\\|[^|]+\\|[^\\s]+ (.+)$"
      for (j in 1:nrow(batch)) {
        if (grepl(sp_pattern, batch$sseqid[j])) {
          hit_descriptions[j] <- gsub(sp_pattern, "\\1", batch$sseqid[j])
        }
      }

      # Extract accession from sseqid
      hit_accessions <- batch$sseqid
      sp_acc_pattern <- "^sp\\|([^|]+)\\|.+$"
      for (j in 1:nrow(batch)) {
        if (grepl(sp_acc_pattern, batch$sseqid[j])) {
          hit_accessions[j] <- gsub(sp_acc_pattern, "\\1", batch$sseqid[j])
        }
      }

      # Prepare statement
      params <- list()
      for (j in 1:nrow(batch)) {
        params <- c(params, list(
          blast_param_id,
          batch$flanking_id[j],
          hit_accessions[j],
          hit_descriptions[j],
          batch$pident[j],
          batch$length[j],
          batch$mismatch[j],
          batch$gapopen[j],
          batch$qstart[j],
          batch$qend[j],
          batch$sstart[j],
          batch$send[j],
          batch$evalue[j],
          batch$bitscore[j]
        ))
      }

      # Execute batch insert
      placeholders <- paste0("(", paste(rep("?", 14), collapse = ", "), ")")
      query <- paste0(
        "INSERT INTO blast_results (blast_param_id, flanking_id, hit_accession, hit_description, ",
        "percent_identity, alignment_length, mismatches, gap_openings, ",
        "query_start, query_end, subject_start, subject_end, e_value, bit_score) VALUES ",
        paste(rep(placeholders, nrow(batch)), collapse = ", ")
      )

      DBI::dbExecute(con, query, params = params)

      # Update progress bar
      if (verbose) {
        pb$update(end_idx / nrow(results))
      }
    }

    # Commit transaction
    DBI::dbExecute(con, "COMMIT")

    if (verbose) message("Imported ", nrow(results), " BLAST results.")

    return(nrow(results))
  }, error = function(e) {
    # Rollback transaction on error
    DBI::dbExecute(con, "ROLLBACK")
    stop("Error importing BLAST results: ", e$message)
  })
}

#' Perform BLAST search on sequences from the database
#'
#' This function retrieves flanking sequences from the database, performs a BLAST search,
#' stores the results back in the database, and captures database metadata for reproducibility.
#'
#' @param con A database connection object.
#' @param project_id The ID of the project to associate with the BLAST search.
#' @param vcf_file_id The ID of the input file containing the VCF data.
#' @param db_path Character string specifying the path to the BLAST database.
#' @param db_name Character string specifying the name of the BLAST database.
#' @param blast_type Character string specifying the type of BLAST search.
#'   Must be either "blastn" or "blastx".
#' @param e_value Numeric value specifying the E-value threshold for hits.
#'   Default is 1e-5.
#' @param max_hits Integer specifying the maximum number of hits to return
#'   per query. Default is 5.
#' @param threads Integer specifying the number of CPU threads to use.
#'   Default is 1.
#' @param output_dir Character string specifying the directory for output files.
#'   Default is the current working directory.
#' @param taxids Optional. Character string specifying NCBI taxonomy IDs to limit the search.
#'   Default is NULL.
#' @param extract_db_metadata Logical. If TRUE, extract and store database metadata.
#'   Default is TRUE.
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE.
#'
#' @return A list containing:
#'   \item{blast_param_id}{The ID of the registered BLAST parameters.}
#'   \item{result_count}{The number of BLAST results imported.}
#'   \item{output_base}{The base name for output files.}
#'   \item{metadata_id}{The ID of the stored database metadata (if extracted).}
#'   \item{db_metadata}{The extracted database metadata (if extracted).}
#'
#' @importFrom DBI dbExecute dbGetQuery
#' @importFrom Biostrings writeXStringSet
#' @export
perform_blast_db <- function(con, project_id, vcf_file_id, db_path, db_name,
                             blast_type = c("blastn", "blastx"),
                             e_value = 1e-5, max_hits = 5, threads = 1,
                             output_dir = getwd(), taxids = NULL,
                             extract_db_metadata = TRUE, verbose = TRUE) {
  # Match arguments
  blast_type <- match.arg(blast_type)

  # Generate output base name
  output_base <- file.path(output_dir, paste0(
    "project_", project_id, "_vcf_", vcf_file_id, "_",
    db_name, "_", blast_type, "_", format(Sys.time(), "%Y%m%d_%H%M%S")
  ))

  # Register BLAST parameters
  if (verbose) message("Registering BLAST parameters...")
  blast_param_id <- register_blast_params(
    con, project_id, blast_type, db_name, db_path, e_value, max_hits, verbose = verbose
  )

  # Extract database metadata before running BLAST
  db_metadata <- NULL
  metadata_id <- NULL

  if (extract_db_metadata) {
    if (verbose) message("Extracting BLAST database metadata...")
    db_metadata <- extract_blast_db_metadata(db_path, db_name, verbose = verbose)

    if (!is.null(db_metadata)) {
      metadata_id <- store_blast_db_metadata(con, blast_param_id, db_metadata, verbose = FALSE)
    }
  }

  # Get flanking sequences
  if (verbose) message("Retrieving flanking sequences...")
  flanking_seqs <- get_flanking_sequences(con, vcf_file_id, as_dna_string_set = TRUE)

  if (length(flanking_seqs) == 0) {
    stop("No flanking sequences found for VCF file ID ", vcf_file_id)
  }

  # Get flanking IDs for mapping results later
  flanking_ids <- DBI::dbGetQuery(
    con,
    "SELECT f.flanking_id, v.chromosome, v.position
     FROM flanking_sequences f
     JOIN vcf_data v ON f.vcf_id = v.vcf_id
     WHERE v.file_id = ?
     ORDER BY v.chromosome, v.position",
    params = list(vcf_file_id)
  )

  # Construct full database path
  db_file <- file.path(db_path, db_name)

  # Check for .nal file (indicator of a multi-volume database)
  is_multi_volume <- file.exists(paste0(db_file, ".nal"))

  if (is_multi_volume) {
    if (verbose) cat("Multi-volume database detected.\n")
  } else {
    # For single-volume databases, check files
    required_extensions <- if(blast_type == "blastn") c(".nhr", ".nin", ".nsq") else c(".phr", ".pin", ".psq")
    missing_files <- sapply(required_extensions, function(ext) !file.exists(paste0(db_file, ext)))
    if (any(missing_files)) {
      stop(paste("Missing required database files:",
                 paste(required_extensions[missing_files], collapse = ", "),
                 "\nPlease check the path and database name."))
    }
  }

  # Write sequences to fasta file
  output_fasta <- paste0(output_base, "_query_seqs.fasta")
  Biostrings::writeXStringSet(flanking_seqs, output_fasta)
  if (verbose) cat("Query sequences written to:", output_fasta, "\n")

  # Set output format
  outfmt <- "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

  # Define output_blast
  output_blast <- paste0(output_base, "_", blast_type, ".txt")

  # Construct blast command
  blast_command <- paste(
    blast_type,
    "-db", shQuote(db_file),
    "-query", shQuote(output_fasta),
    "-out", shQuote(output_blast),
    "-outfmt", shQuote(outfmt),
    "-evalue", e_value,
    "-max_target_seqs", max_hits,
    "-num_threads", threads
  )

  # Add taxids if provided
  if (!is.null(taxids)) {
    blast_command <- paste(blast_command, "-taxids", taxids)
  }

  # Run blast
  if (verbose) {
    cat("Executing BLAST command...\n")
    cat("Command:", blast_command, "\n")

    # Display database info if we have it
    if (!is.null(db_metadata)) {
      cat("\n=== Database Information ===\n")
      cat("Database:", db_metadata$db_title %||% "Unknown", "\n")
      cat("Sequences:", format(db_metadata$num_sequences %||% 0, big.mark = ","), "\n")
      cat("Total length:", format(db_metadata$total_length %||% 0, big.mark = ","), "\n")
      cat("Database date:", db_metadata$db_date %||% "Unknown", "\n")
      if (!is.null(db_metadata$db_version)) {
        cat("Version:", db_metadata$db_version, "\n")
      }
      cat("=============================\n\n")
    }
  }

  system_result <- system(blast_command)

  if (system_result != 0) {
    stop(paste("BLAST command failed with exit status:", system_result))
  }

  if (verbose) cat("BLAST results written to:", output_blast, "\n")

  # Import results into database
  if (verbose) message("Importing BLAST results into database...")

  # Check if results file exists and has content
  if (!file.exists(output_blast) || file.size(output_blast) == 0) {
    if (verbose) message("No BLAST hits found.")
    result_count <- 0
  } else {
    result_count <- import_blast_results(
      con, blast_param_id, output_blast, flanking_ids$flanking_id, verbose = verbose
    )
  }

  # Return comprehensive summary
  result <- list(
    blast_param_id = blast_param_id,
    result_count = result_count,
    output_base = output_base,
    metadata_id = metadata_id,
    db_metadata = db_metadata
  )

  if (verbose) {
    message("\n=== BLAST Search Summary ===")
    message("BLAST parameter ID: ", blast_param_id)
    message("Results imported: ", result_count)
    message("Output base: ", output_base)
    if (!is.null(metadata_id)) {
      message("Database metadata ID: ", metadata_id)
    }
    message("===========================")
  }

  # Update analysis report if it exists
  tryCatch({
    blast_message <- paste0(
      "BLAST search completed: ", result_count, " results from ", db_name, " database",
      if (!is.null(db_metadata)) paste0(" (", format(db_metadata$num_sequences %||% 0, big.mark = ","), " sequences)")
    )

    update_analysis_report(
      con, project_id,
      section = "blast_search",
      message = blast_message,
      verbose = FALSE
    )
  }, error = function(e) {
    # Silently ignore if no report exists
  })

  return(result)
}

#' Get BLAST results from the database
#'
#' This function retrieves BLAST results from the database.
#'
#' @param con A database connection object.
#' @param blast_param_id The ID of the BLAST parameters.
#' @param e_value_threshold Optional. E-value threshold for filtering results. Default is NULL.
#' @param max_hits_per_query Optional. Maximum number of hits to return per query. Default is NULL.
#' @param include_vcf_info Logical. If TRUE, include VCF information. Default is TRUE.
#'
#' @return A data frame containing the BLAST results.
#'
#' @importFrom DBI dbGetQuery
#' @export
get_blast_results <- function(con, blast_param_id, e_value_threshold = NULL,
                          max_hits_per_query = NULL, include_vcf_info = TRUE) {
  # Check if BLAST parameters exist
  params <- DBI::dbGetQuery(
    con,
    "SELECT * FROM blast_parameters WHERE blast_param_id = ?",
    params = list(blast_param_id)
  )

  if (nrow(params) == 0) {
    stop("BLAST parameters with ID ", blast_param_id, " not found.")
  }

  # Build query
  if (include_vcf_info) {
    base_query <- "
      SELECT br.*, f.vcf_id, v.chromosome, v.position, v.ref, v.alt
      FROM blast_results br
      JOIN flanking_sequences f ON br.flanking_id = f.flanking_id
      JOIN vcf_data v ON f.vcf_id = v.vcf_id
      WHERE br.blast_param_id = ?
    "
  } else {
    base_query <- "
      SELECT br.*
      FROM blast_results br
      WHERE br.blast_param_id = ?
    "
  }

  # Add filters
  params_list <- list(blast_param_id)

  if (!is.null(e_value_threshold)) {
    base_query <- paste0(base_query, " AND br.e_value <= ?")
    params_list <- c(params_list, list(e_value_threshold))
  }

  # Order by flanking ID and e-value
  base_query <- paste0(base_query, " ORDER BY br.flanking_id, br.e_value")

  # Execute query
  results <- DBI::dbGetQuery(con, base_query, params = params_list)

  if (nrow(results) == 0) {
    return(data.frame())
  }

  # Apply max hits per query filter if needed
  if (!is.null(max_hits_per_query) && max_hits_per_query > 0) {
    # Group by flanking_id and keep only top hits
    results <- do.call(rbind, lapply(split(results, results$flanking_id), function(group) {
      group[order(group$e_value), ][1:min(nrow(group), max_hits_per_query), ]
    }))

    # Re-order by flanking_id and e_value
    results <- results[order(results$flanking_id, results$e_value), ]
  }

  return(results)
}

#' Count BLAST results in the database
#'
#' @param con A database connection object.
#' @param blast_param_id Optional. The ID of the BLAST parameters. Default is NULL.
#' @param project_id Optional. The ID of the project. Default is NULL.
#'
#' @return The number of BLAST results.
#'
#' @importFrom DBI dbGetQuery
#' @export
count_blast_results <- function(con, blast_param_id = NULL, project_id = NULL) {
  if (!is.null(blast_param_id)) {
    # Count results for specific BLAST parameters
    DBI::dbGetQuery(
      con,
      "SELECT COUNT(*) AS count FROM blast_results WHERE blast_param_id = ?",
      params = list(blast_param_id)
    )$count
  } else if (!is.null(project_id)) {
    # Count results for a specific project
    DBI::dbGetQuery(
      con,
      "SELECT COUNT(*) AS count FROM blast_results r
       JOIN blast_parameters p ON r.blast_param_id = p.blast_param_id
       WHERE p.project_id = ?",
      params = list(project_id)
    )$count
  } else {
    # Count all results
    DBI::dbGetQuery(con, "SELECT COUNT(*) AS count FROM blast_results")$count
  }
}

#' Delete BLAST results from the database
#'
#' @param con A database connection object.
#' @param blast_param_id The ID of the BLAST parameters.
#' @param confirm Logical. If TRUE, user confirmation will be required. Default is TRUE.
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE.
#'
#' @return Logical. TRUE if the deletion was successful.
#'
#' @importFrom DBI dbExecute dbGetQuery
#' @export
delete_blast_results <- function(con, blast_param_id, confirm = TRUE, verbose = TRUE) {
  # Check if BLAST parameters exist
  params <- DBI::dbGetQuery(
    con,
    "SELECT * FROM blast_parameters WHERE blast_param_id = ?",
    params = list(blast_param_id)
  )

  if (nrow(params) == 0) {
    stop("BLAST parameters with ID ", blast_param_id, " not found.")
  }

  # Count results
  count <- DBI::dbGetQuery(
    con,
    "SELECT COUNT(*) AS count FROM blast_results WHERE blast_param_id = ?",
    params = list(blast_param_id)
  )$count

  if (count == 0) {
    if (verbose) message("No BLAST results found for parameter ID ", blast_param_id)
    return(TRUE)
  }

  # Confirm deletion
  if (confirm) {
    answer <- readline(paste0("Are you sure you want to delete ", count,
                             " BLAST results for parameter ID ", blast_param_id,
                             "? (y/n): "))

    if (tolower(answer) != "y") {
      message("BLAST results deletion cancelled.")
      return(FALSE)
    }
  }

  # Start transaction
  DBI::dbExecute(con, "BEGIN TRANSACTION")

  tryCatch({
    # Delete BLAST results
    deleted <- DBI::dbExecute(
      con,
      "DELETE FROM blast_results WHERE blast_param_id = ?",
      params = list(blast_param_id)
    )

    # Commit transaction
    DBI::dbExecute(con, "COMMIT")

    if (verbose) message("Deleted ", deleted, " BLAST results.")

    return(TRUE)
  }, error = function(e) {
    # Rollback transaction on error
    DBI::dbExecute(con, "ROLLBACK")
    stop("Error deleting BLAST results: ", e$message)
  })
}
