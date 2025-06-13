# EXPORTED

#' BLAST search and results management functions for funseqR
#'
#' These functions manage BLAST parameters and results within the funseqR database.
#'

#' Register BLAST parameters in the database
#'
#' This function stores BLAST search parameters in the database.
#'
#' @param con A database connection object.
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
register_blast_params <- function(con, blast_type, db_name, db_path,
                           e_value, max_hits, verbose = TRUE) {

  # Validate blast_type
  blast_type <- match.arg(blast_type, c("blastn", "blastx", "diamond_blastn", "diamond_blastx"))

  # Register BLAST parameters
  current_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  DBI::dbExecute(
    con,
    "INSERT INTO blast_parameters (blast_type, db_name, db_path, e_value, max_hits, execution_date)
     VALUES (?, ?, ?, ?, ?, ?)",
    params = list(blast_type, db_name, db_path, e_value, max_hits, current_time)
  )

  # Get the ID of the newly registered parameters
  param_id <- DBI::dbGetQuery(
    con,
    "SELECT blast_param_id FROM blast_parameters
     WHERE execution_date = ?
     ORDER BY blast_param_id DESC LIMIT 1",
    params = list(current_time)
  )$blast_param_id[1]

  if (verbose) message("Registered BLAST parameters with ID ", param_id)

  return(param_id)
}

#' List all BLAST parameters in the database
#'
#' @param con A database connection object.
#'
#' @return A data frame containing all BLAST parameters in the database.
#'
#' @importFrom DBI dbGetQuery
#' @export
list_blast_params <- function(con) {
  # Get all parameters
  DBI::dbGetQuery(
    con,
    "SELECT * FROM blast_parameters ORDER BY blast_param_id"
  )
}

#' Import BLAST results into the database
#'
#' This function imports BLAST results from a file into the database with proper
#' mapping between query sequence IDs and flanking sequence IDs.
#'
#' @param con A database connection object.
#' @param blast_param_id The ID of the BLAST parameters.
#' @param results_file Path to the BLAST results file.
#' @param vcf_file_id The VCF file ID for looking up flanking sequences.
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE.
#'
#' @return The number of BLAST results imported.
#'
#' @importFrom DBI dbExecute dbGetQuery
#' @importFrom progress progress_bar
#' @export
import_blast_results <- function(con, blast_param_id, results_file, vcf_file_id, verbose = TRUE) {
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

  if (verbose) message("Found ", nrow(results), " BLAST hits to import")

  # Create mapping from BLAST query IDs to flanking sequence IDs
  if (verbose) message("Creating mapping from query IDs to flanking sequence IDs...")

  # Get all flanking sequences with their chromosome and position info
  flanking_lookup <- DBI::dbGetQuery(
    con,
    "SELECT f.flanking_id, v.chromosome, v.position, f.start_position, f.end_position
     FROM flanking_sequences f
     JOIN vcf_data v ON f.vcf_id = v.vcf_id
     WHERE v.file_id = ?
     ORDER BY v.chromosome, v.position",
    params = list(vcf_file_id)
  )

  if (nrow(flanking_lookup) == 0) {
    stop("No flanking sequences found for VCF file ID ", vcf_file_id)
  }

  # Create expected sequence names in the same format as get_flanking_sequences
  flanking_lookup$expected_qseqid <- paste0(
    flanking_lookup$chromosome, ":",
    flanking_lookup$position, " (",
    flanking_lookup$start_position, "-",
    flanking_lookup$end_position, ")"
  )

  # Also create a version without the coordinate suffix for partial matching
  flanking_lookup$partial_qseqid <- paste0(
    flanking_lookup$chromosome, ":",
    flanking_lookup$position
  )

  if (verbose) {
    message("Created lookup table with ", nrow(flanking_lookup), " flanking sequences")
  }

  # Map query IDs to flanking IDs
  unique_queries <- unique(results$qseqid)
  query_to_flanking <- data.frame(
    qseqid = character(0),
    flanking_id = integer(0),
    stringsAsFactors = FALSE
  )

  unmatched_queries <- character(0)

  for (qid in unique_queries) {
    # Try exact match first
    exact_match <- flanking_lookup[flanking_lookup$expected_qseqid == qid, ]

    if (nrow(exact_match) > 0) {
      query_to_flanking <- rbind(query_to_flanking, data.frame(
        qseqid = qid,
        flanking_id = exact_match$flanking_id[1],
        stringsAsFactors = FALSE
      ))
    } else {
      # Try partial match (chromosome:position)
      partial_match <- flanking_lookup[flanking_lookup$partial_qseqid == qid, ]

      if (nrow(partial_match) > 0) {
        query_to_flanking <- rbind(query_to_flanking, data.frame(
          qseqid = qid,
          flanking_id = partial_match$flanking_id[1],
          stringsAsFactors = FALSE
        ))
      } else {
        # Try to extract chromosome and position from qseqid manually
        # Handle formats like "CauratusV1_scaffold_1004:51240" or "LG1:12345"
        parts <- strsplit(qid, ":")[[1]]
        if (length(parts) >= 2) {
          chrom <- parts[1]
          # Extract position (remove any trailing content)
          pos_part <- parts[2]
          pos_match <- regexpr("^[0-9]+", pos_part)
          if (pos_match > 0) {
            position <- as.integer(regmatches(pos_part, pos_match))

            manual_match <- flanking_lookup[
              flanking_lookup$chromosome == chrom &
              flanking_lookup$position == position,
            ]

            if (nrow(manual_match) > 0) {
              query_to_flanking <- rbind(query_to_flanking, data.frame(
                qseqid = qid,
                flanking_id = manual_match$flanking_id[1],
                stringsAsFactors = FALSE
              ))
            } else {
              unmatched_queries <- c(unmatched_queries, qid)
            }
          } else {
            unmatched_queries <- c(unmatched_queries, qid)
          }
        } else {
          unmatched_queries <- c(unmatched_queries, qid)
        }
      }
    }
  }

  if (verbose) {
    message("Successfully mapped ", nrow(query_to_flanking), " out of ", length(unique_queries), " unique query sequences")
    if (length(unmatched_queries) > 0) {
      message("Warning: Could not map ", length(unmatched_queries), " query sequences")
    }
  }

  # Filter results to only include those we can map
  mapped_results <- merge(results, query_to_flanking, by = "qseqid")

  if (nrow(mapped_results) == 0) {
    warning("No BLAST results could be mapped to flanking sequences")
    return(0)
  }

  if (verbose) message("Importing ", nrow(mapped_results), " mapped BLAST results...")

  # Start transaction
  DBI::dbExecute(con, "BEGIN TRANSACTION")

  # Set up progress bar
  if (verbose) {
    pb <- progress::progress_bar$new(
      format = "[:bar] :percent ETA: :eta",
      total = nrow(mapped_results),
      clear = FALSE,
      width = 60
    )
  }

  # Import results in batches
  batch_size <- 1000
  num_batches <- ceiling(nrow(mapped_results) / batch_size)

  tryCatch({
    for (i in 1:num_batches) {
      start_idx <- (i - 1) * batch_size + 1
      end_idx <- min(i * batch_size, nrow(mapped_results))
      batch <- mapped_results[start_idx:end_idx, ]

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
        pb$update(end_idx / nrow(mapped_results))
      }
    }

    # Commit transaction
    DBI::dbExecute(con, "COMMIT")

    if (verbose) message("Successfully imported ", nrow(mapped_results), " BLAST results.")

    return(nrow(mapped_results))
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
#' @param vcf_file_id The ID of the input file containing the VCF data.
#' @param db_path Character string specifying the path to the BLAST database.
#' @param db_name Character string specifying the name of the BLAST database.
#' @param blast_type Character string specifying the type of BLAST search.
#'   Must be either "blastn", "blastx", "diamond_blastn", or "diamond_blastx".
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
#' @param seq_type Type of sequence to search: "raw", "orf_nuc", or "orf_aa". Default is "raw".
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
perform_blast_db <- function(con, vcf_file_id, db_path, db_name,
                             blast_type = c("blastn", "blastx", "diamond_blastn", "diamond_blastx"),
                             e_value = 1e-5, max_hits = 5, threads = 1,
                             output_dir = getwd(), taxids = NULL,
                             extract_db_metadata = TRUE, seq_type = "raw", verbose = TRUE) {
  # Match arguments
  blast_type <- match.arg(blast_type)

  # Determine engine and search type
  is_diamond <- grepl("^diamond_", blast_type)
  search_type <- if (is_diamond) gsub("^diamond_", "", blast_type) else blast_type

  # Generate output base name
  output_base <- file.path(output_dir, paste0(
    "vcf_", vcf_file_id, "_",
    db_name, "_", blast_type, "_", format(Sys.time(), "%Y%m%d_%H%M%S")
  ))

  # Register BLAST parameters
  if (verbose) message("Registering BLAST parameters...")
  blast_param_id <- register_blast_params(
    con, blast_type, db_name, db_path, e_value, max_hits, verbose = verbose
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
  if (verbose) message("Retrieving flanking sequences (", seq_type, ")...")
  flanking_seqs <- get_flanking_sequences(con, vcf_file_id, seq_type = seq_type, as_dna_string_set = TRUE)

  if (length(flanking_seqs) == 0) {
    stop("No flanking sequences found for VCF file ID ", vcf_file_id)
  }


  # Construct full database path
  db_file <- file.path(db_path, db_name)

  # Check database files based on search engine
  if (is_diamond) {
    # For DIAMOND, check for .dmnd file
    dmnd_file <- paste0(db_file, ".dmnd")
    if (!file.exists(dmnd_file)) {
      stop("DIAMOND database file not found: ", dmnd_file,
           "\nPlease create a DIAMOND database or use BLAST instead.")
    }
    if (verbose) cat("DIAMOND database found:", dmnd_file, "\n")
  } else {
    # Check for .nal file (indicator of a multi-volume database)
    is_multi_volume <- file.exists(paste0(db_file, ".nal"))

    if (is_multi_volume) {
      if (verbose) cat("Multi-volume database detected.\n")
    } else {
      # For single-volume databases, check files
      required_extensions <- if(search_type == "blastn") c(".nhr", ".nin", ".nsq") else c(".phr", ".pin", ".psq")
      missing_files <- sapply(required_extensions, function(ext) !file.exists(paste0(db_file, ext)))
      if (any(missing_files)) {
        stop(paste("Missing required database files:",
                   paste(required_extensions[missing_files], collapse = ", "),
                   "\nPlease check the path and database name."))
      }
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

  # Construct command based on search engine
  if (is_diamond) {
    # Construct DIAMOND command
    db_param <- paste0(db_file, ".dmnd")
    blast_command <- paste(
      "diamond", search_type,
      "--db", shQuote(db_param),
      "--query", shQuote(output_fasta),
      "--out", shQuote(output_blast),
      "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
      "qstart", "qend", "sstart", "send", "evalue", "bitscore",
      "--evalue", e_value,
      "--max-target-seqs", max_hits,
      "--threads", threads
    )

    # Add taxids if provided (DIAMOND uses different syntax)
    if (!is.null(taxids)) {
      blast_command <- paste(blast_command, "--taxonlist", taxids)
    }
  } else {
    # Construct BLAST command
    blast_command <- paste(
      search_type,
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
  }

  # Run search
  if (verbose) {
    if (is_diamond) {
      cat("Executing DIAMOND command...\n")
    } else {
      cat("Executing BLAST command...\n")
    }
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
    engine_name <- if (is_diamond) "DIAMOND" else "BLAST"
    stop(paste(engine_name, "command failed with exit status:", system_result))
  }

  engine_name <- if (is_diamond) "DIAMOND" else "BLAST"
  if (verbose) cat(engine_name, "results written to:", output_blast, "\n")

  # Import results into database
  if (verbose) message("Importing ", engine_name, " results into database...")

  # Check if results file exists and has content
  if (!file.exists(output_blast) || file.size(output_blast) == 0) {
    if (verbose) message("No ", engine_name, " hits found.")
    result_count <- 0
  } else {
    result_count <- import_blast_results(
      con, blast_param_id, output_blast, vcf_file_id, verbose = verbose
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
    # Get query sequence count
    query_count <- length(flanking_seqs)
    hit_rate <- if (query_count > 0) round((result_count / query_count) * 100, 1) else 0

    blast_message <- paste0(
      "**BLAST Search Completed**\n\n",
      "- **Search type:** ", blast_type, "\n",
      "- **Database:** ", db_name, "\n",
      if (!is.null(db_metadata)) {
        paste0(
          "- **Database info:** ", db_metadata$db_title %||% "Unknown", "\n",
          "- **Database sequences:** ", format(db_metadata$num_sequences %||% 0, big.mark = ","), "\n",
          "- **Database date:** ", db_metadata$db_date %||% "Unknown", "\n"
        )
      } else {
        "- **Database info:** Not available\n"
      },
      "- **Query sequences:** ", format(query_count, big.mark = ","), "\n",
      "- **BLAST hits found:** ", format(result_count, big.mark = ","), "\n",
      "- **Hit rate:** ", hit_rate, "%\n",
      "- **E-value threshold:** ", e_value, "\n",
      "- **Max hits per query:** ", max_hits, "\n",
      "- **Threads used:** ", threads
    )

    # TODO: Update this when report functions are fixed
    # update_analysis_report(
    #   con,
    #   section = "blast_search",
    #   message = blast_message,
    #   verbose = FALSE
    # )
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
#'
#' @return The number of BLAST results.
#'
#' @importFrom DBI dbGetQuery
#' @export
count_blast_results <- function(con, blast_param_id = NULL) {
  if (!is.null(blast_param_id)) {
    # Count results for specific BLAST parameters
    DBI::dbGetQuery(
      con,
      "SELECT COUNT(*) AS count FROM blast_results WHERE blast_param_id = ?",
      params = list(blast_param_id)
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


# INTERNAL

#' Get BLAST parameters from the database
#'
#' @param con A database connection object.
#' @param blast_param_id The ID of the BLAST parameters to retrieve.
#'
#' @return A data frame containing the BLAST parameters.
#'
#' @importFrom DBI dbGetQuery
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

