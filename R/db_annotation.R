#' Ensure the uniprot_cache table exists
#'
#' @param con A database connection object
#' @param verbose Logical. If TRUE, print progress information
#' @return Invisible NULL
#' @noRd
ensure_uniprot_cache_table <- function(con, verbose = FALSE) {
  # Check if table exists
  tables <- DBI::dbListTables(con)
  if ("uniprot_cache" %in% tables) {
    if (verbose) message("UniProt cache table already exists")
    return(invisible(NULL))
  }

  # Create table
  if (verbose) message("Creating UniProt cache table...")

  DBI::dbExecute(con, "
    CREATE TABLE uniprot_cache (
      cache_id INTEGER PRIMARY KEY,
      accession TEXT NOT NULL UNIQUE,
      response_json TEXT NOT NULL,
      retrieval_date TEXT NOT NULL
    )
  ")

  DBI::dbExecute(con, "CREATE INDEX idx_uniprot_cache_accession ON uniprot_cache (accession)")

  if (verbose) message("UniProt cache table created successfully")
  return(invisible(NULL))
}

#' Get cached UniProt data
#'
#' @param con A database connection object
#' @param accession The UniProt accession number
#' @return The cached UniProt data or NULL if not found
#' @noRd
get_cached_uniprot_data <- function(con, accession) {
  # Check if data exists in cache
  cache_query <- "SELECT response_json FROM uniprot_cache WHERE accession = ?"
  cache_result <- DBI::dbGetQuery(con, cache_query, params = list(accession))

  if (nrow(cache_result) == 0) {
    return(NULL)
  }

  # Parse the JSON response
  json_data <- cache_result$response_json[1]
  data <- jsonlite::fromJSON(json_data, flatten = TRUE)

  # Create a response object
  list(
    accession = accession,
    url = NA_character_,  # We don't store the URL
    status_code = 200,    # Pretend it was a successful request
    content = json_data,
    data = data,
    error = NA_character_,
    from_cache = TRUE
  )
}

#' Store UniProt data in cache
#'
#' @param con A database connection object
#' @param accession The UniProt accession number
#' @param response The UniProt API response
#' @return Invisible NULL
#' @noRd
store_uniprot_data <- function(con, accession, response) {
  # Skip if there was an error
  if (!is.na(response$error)) {
    return(invisible(NULL))
  }

  # Convert data to JSON
  json_data <- response$content

  # Check if this accession already exists
  check_query <- "SELECT cache_id FROM uniprot_cache WHERE accession = ?"
  check_result <- DBI::dbGetQuery(con, check_query, params = list(accession))

  current_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  if (nrow(check_result) > 0) {
    # Update existing record
    update_query <- "UPDATE uniprot_cache SET response_json = ?, retrieval_date = ? WHERE accession = ?"
    DBI::dbExecute(con, update_query, params = list(json_data, current_time, accession))
  } else {
    # Insert new record
    insert_query <- "INSERT INTO uniprot_cache (accession, response_json, retrieval_date) VALUES (?, ?, ?)"
    DBI::dbExecute(con, insert_query, params = list(accession, json_data, current_time))
  }

  return(invisible(NULL))
}

#' Test connection to the UniProt API using working endpoints
#'
#' @param verbose Logical. If TRUE, print progress information
#' @return Logical. TRUE if the connection is successful
#' @export
test_uniprot_connection <- function(verbose = FALSE) {
  # Use a different test approach that matches our successful configurations
  if (verbose) message("Testing UniProt API connection with working endpoints...")

  # Use API client headers based on our tests
  headers <- c(
    "User-Agent" = "funseqR/R API client",
    "Accept" = "application/json"
  )

  # Skip base URL test that we know fails with 403
  # Test the search endpoint directly
  test_acc <- "P99999"  # Cytochrome C - a known protein
  test_url <- paste0("https://rest.uniprot.org/uniprotkb/search?query=accession:",
                     test_acc, "&format=json")

  tryCatch({
    if (verbose) message("Testing REST API search endpoint...")
    response <- httr::GET(test_url, httr::add_headers(.headers = headers))
    status <- httr::status_code(response)

    if (status == 200) {
      if (verbose) message("REST API connection successful!")
      return(TRUE)
    } else {
      if (verbose) message("REST API returned status: ", status, ", trying legacy API...")

      # Try legacy endpoint as fallback
      legacy_url <- paste0("https://www.uniprot.org/uniprot?query=accession:",
                           test_acc, "&format=json")
      legacy_response <- httr::GET(legacy_url, httr::add_headers(.headers = headers))
      legacy_status <- httr::status_code(legacy_response)

      if (legacy_status == 200) {
        if (verbose) message("Legacy API connection successful!")
        return(TRUE)
      } else {
        if (verbose) message("Legacy API returned status: ", legacy_status)
        return(FALSE)
      }
    }
  }, error = function(e) {
    if (verbose) message("Error testing UniProt connection: ", e$message)
    return(FALSE)
  })
}

#' Updated connection check for annotate_blast_results function
#'
#' Place this code within your annotate_blast_results function
#' to replace the existing connection check
check_uniprot_connection <- function(verbose = TRUE) {
  if (verbose) message("Testing UniProt API connection...")
  connection_result <- test_uniprot_connection(verbose = verbose)

  if (!connection_result) {
    message("\nUniProt API connection failed. Possible causes:")
    message("1. Temporary service disruption")
    message("2. Network connectivity issues")
    message("3. API endpoint changes")
    message("\nOptions:")
    message("- Try again later")
    message("- Use offline_mode=TRUE to continue without annotations")
    message("- Check for package updates")
    stop("Cannot connect to UniProt API")
  }

  if (verbose) message("UniProt API connection successful.")
  return(TRUE)
}

#' Extract UniProt information from API response with robust error handling
#'
#' @param uniprot_data The response from query_uniprot_api
#' @param debug If TRUE, print debugging information about the structure of the response
#'
#' @return A list containing extracted information
#'
#' @noRd
extract_uniprot_info <- function(uniprot_data, debug = FALSE) {
  tryCatch({
    if (is.null(uniprot_data$data) || is.null(uniprot_data$data$results) || length(uniprot_data$data$results) == 0) {
      warning("No data found for accession ", uniprot_data$accession)
      return(NULL)
    }

    # Get the first result
    result <- uniprot_data$data$results[[1]]

    # Debug: Print structure if requested
    if (debug) {
      message("Debug info for accession ", uniprot_data$accession, ":")
      message("Available fields: ", paste(names(result), collapse=", "))

      # Check for genes field
      if (!is.null(result$genes)) {
        message("Genes field present with length: ", length(result$genes))
        if (length(result$genes) > 0) {
          message("First gene structure: ")
          print(str(result$genes[[1]]))
        }
      } else {
        message("No genes field found")
      }

      # Check for cross references
      if (!is.null(result$uniProtKBCrossReferences)) {
        message("Cross references present with length: ", length(result$uniProtKBCrossReferences))

        # Look for GO references
        go_refs <- which(sapply(result$uniProtKBCrossReferences, function(r) {
          !is.null(r$database) && r$database == "GO"
        }))

        message("Found ", length(go_refs), " GO references")
        if (length(go_refs) > 0) {
          message("First GO reference structure:")
          print(str(result$uniProtKBCrossReferences[[go_refs[1]]]))
        }

        # Look for KEGG references
        kegg_refs <- which(sapply(result$uniProtKBCrossReferences, function(r) {
          !is.null(r$database) && r$database == "KEGG"
        }))

        message("Found ", length(kegg_refs), " KEGG references")
        if (length(kegg_refs) > 0) {
          message("First KEGG reference structure:")
          print(str(result$uniProtKBCrossReferences[[kegg_refs[1]]]))
        }
      } else {
        message("No uniProtKBCrossReferences field found")
      }
    }

    # Check if result is a proper list structure
    if (is.null(result) || is.atomic(result)) {
      warning("Unexpected result format for accession ", uniprot_data$accession)
      return(create_empty_annotation(uniprot_data$accession))
    }

    # Extract gene names
    gene_names <- ""
    if (!is.null(result$genes) && is.list(result$genes) && length(result$genes) > 0) {
      gene_name_parts <- c()

      for (gene in result$genes) {
        # Add gene name if available
        if (is.list(gene) && !is.null(gene$geneName) && is.list(gene$geneName) && !is.null(gene$geneName$value)) {
          gene_name_parts <- c(gene_name_parts, gene$geneName$value)
        }

        # Add synonyms if available
        if (is.list(gene) && !is.null(gene$synonyms) && is.list(gene$synonyms) && length(gene$synonyms) > 0) {
          for (syn in gene$synonyms) {
            if (is.list(syn) && !is.null(syn$value)) {
              gene_name_parts <- c(gene_name_parts, syn$value)
            }
          }
        }

        # Add ORF names if available
        if (is.list(gene) && !is.null(gene$orfNames) && is.list(gene$orfNames) && length(gene$orfNames) > 0) {
          for (orf in gene$orfNames) {
            if (is.list(orf) && !is.null(orf$value)) {
              gene_name_parts <- c(gene_name_parts, orf$value)
            }
          }
        }
      }

      gene_names <- paste(unique(gene_name_parts), collapse = ";")
    }

    # Extract GO terms
    go_terms <- data.frame(
      go_id = character(0),
      go_term = character(0),
      go_category = character(0),
      go_evidence = character(0),
      stringsAsFactors = FALSE
    )

    # Handle GO terms carefully with explicit list creation
    if (!is.null(result$uniProtKBCrossReferences) && is.list(result$uniProtKBCrossReferences) && length(result$uniProtKBCrossReferences) > 0) {
      go_refs <- list()

      # Find GO references
      for (i in seq_along(result$uniProtKBCrossReferences)) {
        ref <- result$uniProtKBCrossReferences[[i]]
        if (is.list(ref) && !is.null(ref$database) && ref$database == "GO") {
          go_refs <- c(go_refs, list(ref))
        }
      }

      # Process each GO reference
      if (length(go_refs) > 0) {
        for (ref in go_refs) {
          if (!is.null(ref$id) && !is.null(ref$properties) && is.list(ref$properties) && length(ref$properties) > 0) {
            # Extract properties safely
            go_term <- NA_character_
            go_evidence <- NA_character_

            for (prop in ref$properties) {
              if (is.list(prop) && !is.null(prop$key) && !is.null(prop$value)) {
                if (prop$key == "GoTerm") go_term <- prop$value
                if (prop$key == "GoEvidenceType") go_evidence <- prop$value
              }
            }

            # Only add if we found a GO term
            if (!is.na(go_term)) {
              new_row <- data.frame(
                go_id = ref$id,
                go_term = go_term,
                go_category = substr(go_term, 1, 1),
                go_evidence = go_evidence,
                stringsAsFactors = FALSE
              )

              # Append to existing go_terms
              go_terms <- rbind(go_terms, new_row)
            }
          }
        }
      }
    }

    # Extract KEGG references
    kegg_refs <- data.frame(
      kegg_id = character(0),
      pathway_name = character(0),
      stringsAsFactors = FALSE
    )

    # Handle KEGG refs carefully
    if (!is.null(result$uniProtKBCrossReferences) && is.list(result$uniProtKBCrossReferences) && length(result$uniProtKBCrossReferences) > 0) {
      for (i in seq_along(result$uniProtKBCrossReferences)) {
        ref <- result$uniProtKBCrossReferences[[i]]
        if (is.list(ref) && !is.null(ref$database) && ref$database == "KEGG" && !is.null(ref$id)) {
          new_row <- data.frame(
            kegg_id = ref$id,
            pathway_name = NA_character_,
            stringsAsFactors = FALSE
          )
          kegg_refs <- rbind(kegg_refs, new_row)
        }
      }
    }

    # Get primary accession and entry name
    accession <- if (!is.null(result$primaryAccession)) result$primaryAccession else uniprot_data$accession
    entry_name <- if (!is.null(result$uniProtkbId)) result$uniProtkbId else NULL

    # Return the extracted information
    list(
      accession = accession,
      entry_name = entry_name,
      gene_names = gene_names,
      go_terms = go_terms,
      kegg_refs = kegg_refs
    )
  }, error = function(e) {
    warning(paste("Error in extract_uniprot_info for", uniprot_data$accession, ":", e$message))
    return(create_empty_annotation(uniprot_data$accession))
  })
}

#' Create an empty annotation structure
#'
#' @param accession The accession number
#' @return A list with empty annotation structures
#' @noRd
create_empty_annotation <- function(accession) {
  list(
    accession = accession,
    entry_name = NA_character_,
    gene_names = "",
    go_terms = data.frame(
      go_id = character(0),
      go_term = character(0),
      go_category = character(0),
      go_evidence = character(0),
      stringsAsFactors = FALSE
    ),
    kegg_refs = data.frame(
      kegg_id = character(0),
      pathway_name = character(0),
      stringsAsFactors = FALSE
    )
  )
}

#' Store annotation in the database
#'
#' This function stores UniProt annotation in the database.
#'
#' @param con A database connection object.
#' @param blast_result_id The ID of the BLAST result to annotate.
#' @param uniprot_info A list containing UniProt information as returned by extract_uniprot_info.
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE.
#'
#' @return The ID of the newly created annotation.
#'
#' @importFrom DBI dbExecute dbGetQuery
#' @noRd
store_annotation <- function(con, blast_result_id, uniprot_info, verbose = TRUE) {
  if (is.null(uniprot_info)) {
    return(NULL)
  }

  # Check if annotation already exists
  existing <- DBI::dbGetQuery(
    con,
    "SELECT annotation_id FROM annotations
     WHERE blast_result_id = ? AND uniprot_accession = ?",
    params = list(blast_result_id, uniprot_info$accession)
  )

  if (nrow(existing) > 0) {
    if (verbose) message("Annotation already exists with ID ", existing$annotation_id[1])
    return(existing$annotation_id[1])
  }

  # Add annotation
  current_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  DBI::dbExecute(
    con,
    "INSERT INTO annotations (blast_result_id, uniprot_accession, entry_name, gene_names, retrieval_date)
     VALUES (?, ?, ?, ?, ?)",
    params = list(
      blast_result_id,
      uniprot_info$accession,
      uniprot_info$entry_name,
      uniprot_info$gene_names,
      current_time
    )
  )

  # Get the ID of the newly created annotation
  annotation_id <- DBI::dbGetQuery(
    con,
    "SELECT annotation_id FROM annotations
     WHERE blast_result_id = ? AND uniprot_accession = ?",
    params = list(blast_result_id, uniprot_info$accession)
  )$annotation_id[1]

  # Add GO terms
  if (nrow(uniprot_info$go_terms) > 0) {
    for (i in 1:nrow(uniprot_info$go_terms)) {
      DBI::dbExecute(
        con,
        "INSERT INTO go_terms (annotation_id, go_id, go_term, go_category, go_evidence)
         VALUES (?, ?, ?, ?, ?)",
        params = list(
          annotation_id,
          uniprot_info$go_terms$go_id[i],
          uniprot_info$go_terms$go_term[i],
          uniprot_info$go_terms$go_category[i],
          uniprot_info$go_terms$go_evidence[i]
        )
      )
    }
  }

  # Add KEGG references
  if (nrow(uniprot_info$kegg_refs) > 0) {
    for (i in 1:nrow(uniprot_info$kegg_refs)) {
      DBI::dbExecute(
        con,
        "INSERT INTO kegg_references (annotation_id, kegg_id, pathway_name)
         VALUES (?, ?, ?)",
        params = list(
          annotation_id,
          uniprot_info$kegg_refs$kegg_id[i],
          uniprot_info$kegg_refs$pathway_name[i]
        )
      )
    }
  }

  return(annotation_id)
}

#' Process BLAST results and retrieve UniProt annotations
#'
#' This function processes BLAST results from the database and retrieves
#' UniProt annotations for the hit accessions.
#'
#' @param con A database connection object.
#' @param blast_param_id The ID of the BLAST parameters.
#' @param max_hits Maximum number of hits to process per query. Default is 5.
#' @param e_value_threshold E-value threshold for filtering BLAST hits. Default is 1e-10.
#' @param batch_size Number of annotations to process in each transaction. Default is 500.
#' @param delay Delay between API calls in seconds. Default is 1.
#' @param offline_mode If TRUE, skips UniProt API and uses basic annotations. Default is FALSE.
#' @param use_cache If TRUE, uses cached API responses if available. Default is FALSE.
#' @param store_cache If TRUE, stores API responses in the database. Default is FALSE.
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE.
#' @param debug_accessions Vector of accessions to debug. Default is NULL.
#'
#' @return A list containing annotation statistics.
#'
#' @importFrom DBI dbExecute dbGetQuery dbListTables
#' @importFrom httr GET add_headers status_code content
#' @importFrom jsonlite fromJSON
#' @export
annotate_blast_results <- function(con, blast_param_id, max_hits = 5, e_value_threshold = 1e-10,
                                   batch_size = 500, delay = 1, offline_mode = FALSE,
                                   use_cache = FALSE, store_cache = FALSE, verbose = TRUE,
                                   debug_accessions = NULL) {
  # Check connection validity
  if (!DBI::dbIsValid(con)) {
    stop("Database connection is invalid. Please reconnect to the database.")
  }

  # Check if BLAST parameters exist
  params <- DBI::dbGetQuery(
    con,
    "SELECT * FROM blast_parameters WHERE blast_param_id = ?",
    params = list(blast_param_id)
  )

  if (nrow(params) == 0) {
    stop("BLAST parameters with ID ", blast_param_id, " not found.")
  }

  # Create or verify cache table if needed
  if (use_cache || store_cache) {
    ensure_uniprot_cache_table(con, verbose = verbose)
  }

  # Get BLAST results
  if (verbose) message("Retrieving BLAST results...")

  # First get unique hit accessions to minimize API calls
  hit_accessions_query <- "
    SELECT DISTINCT hit_accession
    FROM blast_results
    WHERE blast_param_id = ? AND e_value <= ?
  "

  hit_accessions <- DBI::dbGetQuery(
    con,
    hit_accessions_query,
    params = list(blast_param_id, e_value_threshold)
  )$hit_accession

  if (length(hit_accessions) == 0) {
    stop("No BLAST hits found for parameter ID ", blast_param_id, " with e-value <= ", e_value_threshold)
  }

  if (verbose) message("Found ", length(hit_accessions), " unique hit accessions to annotate.")

  # Check if we should use offline mode
  use_api <- !offline_mode

  # Test UniProt API connection only if not in offline mode and not exclusively using cache
  if (use_api && !(use_cache && !store_cache)) {
    if (verbose) message("Testing UniProt API connection...")
    # Rest of API testing code remains unchanged...
    check_uniprot_connection(verbose = verbose)
  } else if (use_cache && !store_cache) {
    if (verbose) message("Using cache only - skipping API connection test")
  }

  # Process in online or offline mode
  if (use_api) {
    # ONLINE MODE: Query UniProt for each accession or use cache
    if (verbose) message("Retrieving UniProt data for accessions...")

    uniprot_data <- list()
    cache_hits <- 0
    api_calls <- 0

    # Create directory for debug output if needed
    if (!is.null(debug_accessions) && length(debug_accessions) > 0) {
      debug_dir <- "uniprot_debug"
      if (!dir.exists(debug_dir)) dir.create(debug_dir)
      if (verbose) message("Debug mode enabled for accessions: ", paste(debug_accessions, collapse=", "))
    }

    for (i in seq_along(hit_accessions)) {
      acc <- hit_accessions[i]
      cache_data <- NULL

      # Check if this accession should be debugged
      enable_debug <- !is.null(debug_accessions) && acc %in% debug_accessions

      # Check cache first if enabled
      if (use_cache) {
        cache_data <- get_cached_uniprot_data(con, acc)
        if (!is.null(cache_data)) {
          uniprot_data[[acc]] <- cache_data
          cache_hits <- cache_hits + 1

          # Show progress
          if (verbose && (i %% 10 == 0 || i == length(hit_accessions))) {
            message("Processed ", i, " of ", length(hit_accessions), " accessions (",
                    cache_hits, " from cache, ", api_calls, " API calls)")
          }

          next  # Skip API call
        }
      }

      # Query UniProt API with debugging if requested
      result <- query_uniprot_api(acc, debug = enable_debug)
      uniprot_data[[acc]] <- result
      api_calls <- api_calls + 1

      # Store in cache if enabled
      if (store_cache) {
        store_uniprot_data(con, acc, result)
      }

      # Show progress
      if (verbose && (i %% 10 == 0 || i == length(hit_accessions))) {
        message("Processed ", i, " of ", length(hit_accessions), " accessions (",
                cache_hits, " from cache, ", api_calls, " API calls)")
      }

      # Add delay between API calls
      if (api_calls > 0 && i < length(hit_accessions)) {
        Sys.sleep(delay)
      }
    }

    if (verbose) {
      message("UniProt data retrieval complete: ", cache_hits, " from cache, ",
              api_calls, " from API")
    }

    # Extract information from UniProt responses
    if (verbose) message("Extracting information from UniProt responses...")

    uniprot_info <- list()
    for (acc in names(uniprot_data)) {
      # Pass debug flag to extract_uniprot_info for selected accessions
      enable_debug <- !is.null(debug_accessions) && acc %in% debug_accessions
      uniprot_info[[acc]] <- extract_uniprot_info(uniprot_data[[acc]], debug = enable_debug)
    }

    valid_info <- uniprot_info[!sapply(uniprot_info, is.null)]

    if (verbose) message("Successfully extracted information for ", length(valid_info), " of ", length(hit_accessions), " accessions.")
  } else {
    # OFFLINE MODE: Skip API queries and create empty info
    if (verbose) message("Skipping UniProt API queries in offline mode.")
    valid_info <- list()
  }

  # Rest of the function remains unchanged...
  # Get BLAST results to annotate
  blast_results_query <- "
    SELECT blast_result_id, flanking_id, hit_accession, e_value
    FROM blast_results
    WHERE blast_param_id = ? AND e_value <= ?
    ORDER BY flanking_id, e_value
  "

  blast_results <- DBI::dbGetQuery(
    con,
    blast_results_query,
    params = list(blast_param_id, e_value_threshold)
  )

  # Limit to max_hits per flanking sequence
  if (!is.null(max_hits) && max_hits > 0) {
    blast_results <- do.call(rbind, lapply(split(blast_results, blast_results$flanking_id), function(group) {
      group[order(group$e_value), ][1:min(nrow(group), max_hits), ]
    }))
  }

  if (verbose) message("Processing ", nrow(blast_results), " BLAST results in batches of ", batch_size)

  # Helper function for basic annotation when API is unavailable or data is missing
  store_basic_annotation <- function(con, blast_result_id, hit_accession, gene_name = NULL, verbose = FALSE) {
    # ... (existing code)
  }

  # Process BLAST results in batches
  successful_annotations <- 0
  annotation_ids <- c()

  for (batch_start in seq(1, nrow(blast_results), by = batch_size)) {
    batch_end <- min(batch_start + batch_size - 1, nrow(blast_results))
    current_batch <- blast_results[batch_start:batch_end, ]

    if (verbose) message("Processing batch ", ceiling(batch_start/batch_size), " of ", ceiling(nrow(blast_results)/batch_size),
                         " (rows ", batch_start, " to ", batch_end, ")")

    # Start transaction for this batch
    DBI::dbExecute(con, "BEGIN TRANSACTION")

    batch_success <- 0
    batch_ids <- c()

    tryCatch({
      for (i in 1:nrow(current_batch)) {
        blast_result_id <- current_batch$blast_result_id[i]
        hit_accession <- current_batch$hit_accession[i]

        if (use_api && !is.null(valid_info[[hit_accession]])) {
          # ONLINE MODE: Store complete annotation
          annotation_id <- store_annotation(
            con,
            blast_result_id,
            valid_info[[hit_accession]],
            verbose = FALSE
          )
        } else {
          # OFFLINE MODE or missing info: Store basic annotation
          annotation_id <- store_basic_annotation(
            con,
            blast_result_id,
            hit_accession,
            verbose = FALSE
          )
        }

        if (!is.null(annotation_id)) {
          batch_success <- batch_success + 1
          batch_ids <- c(batch_ids, annotation_id)
        }
      }

      # Commit transaction for this batch
      DBI::dbExecute(con, "COMMIT")

      # Update totals
      successful_annotations <- successful_annotations + batch_success
      annotation_ids <- c(annotation_ids, batch_ids)

      if (verbose) message("Batch completed: ", batch_success, " annotations stored")

    }, error = function(e) {
      # Rollback transaction on error
      DBI::dbExecute(con, "ROLLBACK")
      warning("Error in batch processing: ", e$message, ". Batch skipped.")
    })
  }

  if (verbose) message("Successfully annotated ", successful_annotations, " of ", nrow(blast_results), " BLAST results.")

  # Count annotated GO terms and KEGG pathways
  go_count <- 0
  kegg_count <- 0

  if (length(annotation_ids) > 0) {
    # Process in chunks of 1000 to avoid too-long SQL queries with many IDs
    go_count <- 0
    kegg_count <- 0

    for (chunk_start in seq(1, length(annotation_ids), by = 1000)) {
      chunk_end <- min(chunk_start + 999, length(annotation_ids))
      chunk_ids <- annotation_ids[chunk_start:chunk_end]

      go_chunk <- DBI::dbGetQuery(
        con,
        paste0("SELECT COUNT(*) AS count FROM go_terms WHERE annotation_id IN (",
               paste(chunk_ids, collapse = ","), ")")
      )$count

      kegg_chunk <- DBI::dbGetQuery(
        con,
        paste0("SELECT COUNT(*) AS count FROM kegg_references WHERE annotation_id IN (",
               paste(chunk_ids, collapse = ","), ")")
      )$count

      go_count <- go_count + go_chunk
      kegg_count <- kegg_count + kegg_chunk
    }
  }

  if (verbose) {
    message("Stored ", go_count, " GO terms and ", kegg_count, " KEGG references.")

    if (use_cache || store_cache) {
      cache_stats <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS count FROM uniprot_cache")$count
      message("UniProt cache now contains ", cache_stats, " entries")
    }

    if (!use_api) {
      message("\nNOTE: Annotation was performed in offline mode.")
      message("You can re-run the annotation later with UniProt API access to get GO terms and KEGG pathways.")
    }
  }

  # Return summary
  list(
    blast_param_id = blast_param_id,
    unique_accessions = length(hit_accessions),
    successful_extractions = if (use_api) length(valid_info) else 0,
    annotated_results = successful_annotations,
    go_terms = go_count,
    kegg_refs = kegg_count,
    offline_mode = !use_api,
    cache_used = use_cache,
    cache_stored = store_cache
  )
}

#' Query the UniProt API for a protein accession using working endpoints
#'
#' @param accession The UniProt accession number
#' @param fields The fields to retrieve from UniProt
#' @param use_rest_api If TRUE, use REST API, otherwise use legacy API
#' @param debug If TRUE, save the API response to a file for debugging
#' @return A list containing API response
#' @noRd
query_uniprot_api <- function(accession,
                              fields = "accession,id,gene_names,go_id,go,xref_kegg",
                              use_rest_api = TRUE,
                              debug = FALSE) {
  # Use API client headers, which we know work based on our tests
  headers <- c(
    "User-Agent" = paste0("funseqR/R API client"),
    "Accept" = "application/json"
  )

  if (use_rest_api) {
    # Use the REST API search endpoint directly (skip base URL test)
    base_url <- "https://rest.uniprot.org/uniprotkb/search"
    query <- paste0("accession:", accession)
    url <- paste0(base_url, "?query=", URLencode(query),
                  "&fields=", URLencode(fields), "&format=json")

    response <- httr::GET(url, httr::add_headers(.headers = headers))
  } else {
    # Fallback to legacy API which also works
    base_url <- "https://www.uniprot.org/uniprot"
    query <- paste0("accession:", accession)
    url <- paste0(base_url, "?query=", URLencode(query),
                  "&columns=", URLencode(fields), "&format=json")

    response <- httr::GET(url, httr::add_headers(.headers = headers))
  }

  # Process response
  if (httr::status_code(response) == 200) {
    content <- httr::content(response, "text", encoding = "UTF-8")
    data <- jsonlite::fromJSON(content, flatten = TRUE)

    # Debug: Save the response if requested
    if (debug) {
      debug_dir <- "uniprot_debug"
      if (!dir.exists(debug_dir)) dir.create(debug_dir)
      debug_file <- file.path(debug_dir, paste0("uniprot_", accession, ".json"))
      writeLines(content, debug_file)
      message("Saved debug response to ", debug_file)
    }

    return(list(
      accession = accession,
      url = url,
      status_code = httr::status_code(response),
      content = content,
      data = data,
      error = NA_character_
    ))
  } else {
    # Try fallback if REST API failed
    if (use_rest_api) {
      return(query_uniprot_api(accession, fields, use_rest_api = FALSE, debug = debug))
    }

    # Both methods failed
    warning(paste("Failed to retrieve data for", accession,
                  "- Status code:", httr::status_code(response)))

    return(list(
      accession = accession,
      url = url,
      status_code = httr::status_code(response),
      content = NA_character_,
      data = NULL,
      error = paste("Failed to retrieve data - Status code:",
                    httr::status_code(response))
    ))
  }
}

#' Get annotations from the database
#'
#' This function retrieves annotations from the database.
#'
#' @param con A database connection object.
#' @param blast_param_id The ID of the BLAST parameters.
#' @param include_go Logical. If TRUE, include GO terms. Default is TRUE.
#' @param include_kegg Logical. If TRUE, include KEGG references. Default is TRUE.
#' @param include_vcf_info Logical. If TRUE, include VCF information. Default is TRUE.
#'
#' @return A list containing:
#'   \item{annotations}{A data frame of annotations.}
#'   \item{go_terms}{A data frame of GO terms (if include_go is TRUE).}
#'   \item{kegg_refs}{A data frame of KEGG references (if include_kegg is TRUE).}
#'
#' @importFrom DBI dbGetQuery
#' @export
get_annotations <- function(con, blast_param_id, include_go = TRUE, include_kegg = TRUE, include_vcf_info = TRUE) {
  # Check if BLAST parameters exist
  params <- DBI::dbGetQuery(
    con,
    "SELECT * FROM blast_parameters WHERE blast_param_id = ?",
    params = list(blast_param_id)
  )

  if (nrow(params) == 0) {
    stop("BLAST parameters with ID ", blast_param_id, " not found.")
  }

  # Build annotations query
  if (include_vcf_info) {
    annotations_query <- "
      SELECT a.*, br.flanking_id, br.hit_accession, br.e_value, br.bit_score,
             f.vcf_id, v.chromosome, v.position, v.ref, v.alt
      FROM annotations a
      JOIN blast_results br ON a.blast_result_id = br.blast_result_id
      JOIN flanking_sequences f ON br.flanking_id = f.flanking_id
      JOIN vcf_data v ON f.vcf_id = v.vcf_id
      WHERE br.blast_param_id = ?
      ORDER BY v.chromosome, v.position, br.e_value
    "
  } else {
    annotations_query <- "
      SELECT a.*, br.flanking_id, br.hit_accession, br.e_value, br.bit_score
      FROM annotations a
      JOIN blast_results br ON a.blast_result_id = br.blast_result_id
      WHERE br.blast_param_id = ?
      ORDER BY br.flanking_id, br.e_value
    "
  }

  # Get annotations
  annotations <- DBI::dbGetQuery(con, annotations_query, params = list(blast_param_id))

  if (nrow(annotations) == 0) {
    return(list(
      annotations = data.frame(),
      go_terms = data.frame(),
      kegg_refs = data.frame()
    ))
  }

  # Get GO terms if requested
  go_terms <- data.frame()
  if (include_go) {
    go_query <- "
      SELECT g.*, a.uniprot_accession, br.flanking_id
      FROM go_terms g
      JOIN annotations a ON g.annotation_id = a.annotation_id
      JOIN blast_results br ON a.blast_result_id = br.blast_result_id
      WHERE br.blast_param_id = ?
      ORDER BY br.flanking_id, a.uniprot_accession, g.go_category, g.go_id
    "

    go_terms <- DBI::dbGetQuery(con, go_query, params = list(blast_param_id))
  }

  # Get KEGG references if requested
  kegg_refs <- data.frame()
  if (include_kegg) {
    kegg_query <- "
      SELECT k.*, a.uniprot_accession, br.flanking_id
      FROM kegg_references k
      JOIN annotations a ON k.annotation_id = a.annotation_id
      JOIN blast_results br ON a.blast_result_id = br.blast_result_id
      WHERE br.blast_param_id = ?
      ORDER BY br.flanking_id, a.uniprot_accession, k.kegg_id
    "

    kegg_refs <- DBI::dbGetQuery(con, kegg_query, params = list(blast_param_id))
  }

  # Return results
  list(
    annotations = annotations,
    go_terms = go_terms,
    kegg_refs = kegg_refs
  )
}

#' Count annotations in the database
#'
#' @param con A database connection object.
#' @param blast_param_id Optional. The ID of the BLAST parameters. Default is NULL.
#' @param project_id Optional. The ID of the project. Default is NULL.
#'
#' @return A list containing counts of annotations, GO terms, and KEGG references.
#'
#' @importFrom DBI dbGetQuery
#' @export
count_annotations <- function(con, blast_param_id = NULL, project_id = NULL) {
  if (!is.null(blast_param_id)) {
    # Count annotations for specific BLAST parameters
    anno_count <- DBI::dbGetQuery(
      con,
      "SELECT COUNT(*) AS count FROM annotations a
       JOIN blast_results br ON a.blast_result_id = br.blast_result_id
       WHERE br.blast_param_id = ?",
      params = list(blast_param_id)
    )$count

    go_count <- DBI::dbGetQuery(
      con,
      "SELECT COUNT(*) AS count FROM go_terms g
       JOIN annotations a ON g.annotation_id = a.annotation_id
       JOIN blast_results br ON a.blast_result_id = br.blast_result_id
       WHERE br.blast_param_id = ?",
      params = list(blast_param_id)
    )$count

    kegg_count <- DBI::dbGetQuery(
      con,
      "SELECT COUNT(*) AS count FROM kegg_references k
       JOIN annotations a ON k.annotation_id = a.annotation_id
       JOIN blast_results br ON a.blast_result_id = br.blast_result_id
       WHERE br.blast_param_id = ?",
      params = list(blast_param_id)
    )$count
  } else if (!is.null(project_id)) {
    # Count annotations for a specific project
    anno_count <- DBI::dbGetQuery(
      con,
      "SELECT COUNT(*) AS count FROM annotations a
       JOIN blast_results br ON a.blast_result_id = br.blast_result_id
       JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
       WHERE bp.project_id = ?",
      params = list(project_id)
    )$count

    go_count <- DBI::dbGetQuery(
      con,
      "SELECT COUNT(*) AS count FROM go_terms g
       JOIN annotations a ON g.annotation_id = a.annotation_id
       JOIN blast_results br ON a.blast_result_id = br.blast_result_id
       JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
       WHERE bp.project_id = ?",
      params = list(project_id)
    )$count

    kegg_count <- DBI::dbGetQuery(
      con,
      "SELECT COUNT(*) AS count FROM kegg_references k
       JOIN annotations a ON k.annotation_id = a.annotation_id
       JOIN blast_results br ON a.blast_result_id = br.blast_result_id
       JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
       WHERE bp.project_id = ?",
      params = list(project_id)
    )$count
  } else {
    # Count all annotations
    anno_count <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS count FROM annotations")$count
    go_count <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS count FROM go_terms")$count
    kegg_count <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS count FROM kegg_references")$count
  }

  list(
    annotations = anno_count,
    go_terms = go_count,
    kegg_refs = kegg_count
  )
}

#' Delete annotations from the database
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
delete_annotations <- function(con, blast_param_id, confirm = TRUE, verbose = TRUE) {
  # Check if BLAST parameters exist
  params <- DBI::dbGetQuery(
    con,
    "SELECT * FROM blast_parameters WHERE blast_param_id = ?",
    params = list(blast_param_id)
  )

  if (nrow(params) == 0) {
    stop("BLAST parameters with ID ", blast_param_id, " not found.")
  }

  # Count annotations
  counts <- count_annotations(con, blast_param_id)

  if (counts$annotations == 0) {
    if (verbose) message("No annotations found for BLAST parameter ID ", blast_param_id)
    return(TRUE)
  }

  # Confirm deletion
  if (confirm) {
    answer <- readline(paste0("Are you sure you want to delete ", counts$annotations,
                             " annotations, ", counts$go_terms, " GO terms, and ",
                             counts$kegg_refs, " KEGG references for BLAST parameter ID ",
                             blast_param_id, "? (y/n): "))

    if (tolower(answer) != "y") {
      message("Annotations deletion cancelled.")
      return(FALSE)
    }
  }

  # Start transaction
  DBI::dbExecute(con, "BEGIN TRANSACTION")

  tryCatch({
    # Get all annotation IDs
    annotation_ids <- DBI::dbGetQuery(
      con,
      "SELECT a.annotation_id FROM annotations a
       JOIN blast_results br ON a.blast_result_id = br.blast_result_id
       WHERE br.blast_param_id = ?",
      params = list(blast_param_id)
    )$annotation_id

    if (length(annotation_ids) > 0) {
      # Delete GO terms
      if (verbose) message("Deleting GO terms...")

      DBI::dbExecute(
        con,
        paste0("DELETE FROM go_terms WHERE annotation_id IN (",
               paste(annotation_ids, collapse = ","), ")")
      )

      # Delete KEGG references
      if (verbose) message("Deleting KEGG references...")

      DBI::dbExecute(
        con,
        paste0("DELETE FROM kegg_references WHERE annotation_id IN (",
               paste(annotation_ids, collapse = ","), ")")
      )

      # Delete annotations
      if (verbose) message("Deleting annotations...")

      DBI::dbExecute(
        con,
        paste0("DELETE FROM annotations WHERE annotation_id IN (",
               paste(annotation_ids, collapse = ","), ")")
      )
    }

    # Commit transaction
    DBI::dbExecute(con, "COMMIT")

    if (verbose) message("Deleted ", counts$annotations, " annotations, ",
                        counts$go_terms, " GO terms, and ",
                        counts$kegg_refs, " KEGG references.")

    return(TRUE)
  }, error = function(e) {
    # Rollback transaction on error
    DBI::dbExecute(con, "ROLLBACK")
    stop("Error deleting annotations: ", e$message)
  })
}
