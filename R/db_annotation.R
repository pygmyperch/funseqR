#' Ensure the uniprot_cache table exists
#'
#' @param con A database connection object
#' @param verbose Logical. If TRUE, print progress information
#' @return Invisible NULL
#' @export
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
#' @export
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
#' @export
store_uniprot_data <- function(con, accession, response, debug = FALSE) {
  if (debug) {
    message("Attempting to store data for accession: ", accession)
    message("Response structure: ", paste(names(response), collapse=", "))
    if (!is.null(response$status_code)) {
      message("Response status code: ", response$status_code)
    }
  }

  # Check if the database is connected
  if (!DBI::dbIsValid(con)) {
    if (debug) message("Database connection is invalid")
    return(invisible(NULL))
  }

  # Prepare JSON data based on what's available in the response
  json_data <- NULL

  # First try to use content directly if it exists and is a string
  if (!is.null(response$content) && is.character(response$content) && nchar(response$content) > 0) {
    json_data <- response$content
    if (debug) message("Using content as JSON data (length: ", nchar(json_data), ")")
  }
  # Next try using the 'data' field if it exists and is a list
  else if (!is.null(response$data) && is.list(response$data)) {
    tryCatch({
      json_data <- jsonlite::toJSON(response$data, auto_unbox = TRUE)
      if (debug) message("Converted data to JSON (length: ", nchar(json_data), ")")
    }, error = function(e) {
      if (debug) message("Failed to convert data to JSON: ", e$message)
    })
  }
  # If we have a successful status code but no useful content, create a minimal JSON
  else if (!is.null(response$status_code) && response$status_code == 200) {
    # Create a minimal JSON with accession and status
    minimal_data <- list(
      accession = accession,
      status = "success",
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    )
    json_data <- jsonlite::toJSON(minimal_data, auto_unbox = TRUE)
    if (debug) message("Created minimal JSON data (length: ", nchar(json_data), ")")
  }

  # If we still don't have JSON data, create an empty entry
  if (is.null(json_data) || nchar(json_data) == 0) {
    if (debug) message("No valid JSON data available, creating empty entry")
    json_data <- paste0('{"accession":"', accession, '","empty":true}')
  }

  # Store the data regardless of content (we'll at least have the accession)
  tryCatch({
    current_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

    # Check if this accession already exists
    check_query <- "SELECT COUNT(*) AS count FROM uniprot_cache WHERE accession = ?"
    check_result <- DBI::dbGetQuery(con, check_query, params = list(accession))

    if (check_result$count > 0) {
      # Update existing record
      if (debug) message("Updating existing record")
      update_query <- "UPDATE uniprot_cache SET response_json = ?, retrieval_date = ? WHERE accession = ?"
      result <- DBI::dbExecute(con, update_query, params = list(json_data, current_time, accession))
      if (debug) message("Updated ", result, " records")
    } else {
      # Insert new record
      if (debug) message("Inserting new record")
      insert_query <- "INSERT INTO uniprot_cache (accession, response_json, retrieval_date) VALUES (?, ?, ?)"
      result <- DBI::dbExecute(con, insert_query, params = list(accession, json_data, current_time))
      if (debug) message("Inserted ", result, " records")
    }

    # Verify the insertion
    verify_query <- "SELECT COUNT(*) AS count FROM uniprot_cache WHERE accession = ?"
    verify_result <- DBI::dbGetQuery(con, verify_query, params = list(accession))

    if (debug) message("Verification: ", verify_result$count, " records found for ", accession)

    return(invisible(NULL))
  }, error = function(e) {
    if (debug) message("Database error: ", e$message)
    return(invisible(NULL))
  })
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

#' Extract UniProt information from API response
#'
#' @param uniprot_data The response from query_uniprot_api
#' @param debug If TRUE, print debugging information
#'
#' @return A list containing extracted information
#' @export
extract_uniprot_info <- function(uniprot_data, debug = FALSE) {
  # Set up result list with empty data
  accession <- if (!is.null(uniprot_data$accession)) uniprot_data$accession else "Unknown"

  result <- list(
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

  # Exit early if there's an error or no data
  if (!is.null(uniprot_data$error) && !is.na(uniprot_data$error)) {
    if (debug) message("Error in UniProt data: ", uniprot_data$error)
    return(result)
  }

  if (is.null(uniprot_data$data) || !is.list(uniprot_data$data)) {
    if (debug) message("No data available")
    return(result)
  }

  if (is.null(uniprot_data$data$results) || length(uniprot_data$data$results) == 0) {
    if (debug) message("No results found")
    return(result)
  }

  # Get the first result
  entry <- uniprot_data$data$results[[1]]

  # Extract basic info
  if (!is.null(entry$primaryAccession))
    result$accession <- entry$primaryAccession

  if (!is.null(entry$uniProtkbId))
    result$entry_name <- entry$uniProtkbId

  # Extract gene names
  if (!is.null(entry$genes) && length(entry$genes) > 0) {
    gene_names <- c()

    for (gene in entry$genes) {
      if (!is.null(gene$geneName) && !is.null(gene$geneName$value)) {
        gene_names <- c(gene_names, gene$geneName$value)
      }
    }

    if (length(gene_names) > 0) {
      result$gene_names <- paste(gene_names, collapse = ";")
    }
  }

  # Extract GO terms
  if (!is.null(entry$uniProtKBCrossReferences)) {
    # Count GO terms first
    go_count <- 0
    kegg_count <- 0

    for (ref in entry$uniProtKBCrossReferences) {
      if (!is.null(ref$database)) {
        if (ref$database == "GO") go_count <- go_count + 1
        if (ref$database == "KEGG") kegg_count <- kegg_count + 1
      }
    }

    if (debug) {
      message("Found ", go_count, " GO references and ", kegg_count, " KEGG references")
    }

    # Extract GO terms
    if (go_count > 0) {
      go_terms <- data.frame(
        go_id = character(go_count),
        go_term = character(go_count),
        go_category = character(go_count),
        go_evidence = character(go_count),
        stringsAsFactors = FALSE
      )

      go_idx <- 1
      for (ref in entry$uniProtKBCrossReferences) {
        if (!is.null(ref$database) && ref$database == "GO" && !is.null(ref$id)) {
          go_id <- ref$id
          go_term <- NA_character_
          go_evidence <- NA_character_

          if (!is.null(ref$properties)) {
            for (prop in ref$properties) {
              if (!is.null(prop$key) && !is.null(prop$value)) {
                if (prop$key == "GoTerm") go_term <- prop$value
                if (prop$key == "GoEvidenceType") go_evidence <- prop$value
              }
            }
          }

          if (!is.na(go_term)) {
            if (debug) message("Adding GO term: ", go_id, " - ", go_term)
            go_terms[go_idx, ] <- list(
              go_id,
              go_term,
              substr(go_term, 1, 1),
              go_evidence
            )
            go_idx <- go_idx + 1
          }
        }
      }

      # Trim any unused rows
      if (go_idx <= go_count) {
        go_terms <- go_terms[1:(go_idx-1), , drop = FALSE]
      }

      result$go_terms <- go_terms
    }

    # Extract KEGG references
    if (kegg_count > 0) {
      kegg_refs <- data.frame(
        kegg_id = character(kegg_count),
        pathway_name = character(kegg_count),
        stringsAsFactors = FALSE
      )

      kegg_idx <- 1
      for (ref in entry$uniProtKBCrossReferences) {
        if (!is.null(ref$database) && ref$database == "KEGG" && !is.null(ref$id)) {
          kegg_id <- ref$id
          pathway_name <- NA_character_

          if (!is.null(ref$properties)) {
            for (prop in ref$properties) {
              if (!is.null(prop$key) && !is.null(prop$value)) {
                if (prop$key == "Description") pathway_name <- prop$value
              }
            }
          }

          if (debug) message("Adding KEGG reference: ", kegg_id)
          kegg_refs[kegg_idx, ] <- list(kegg_id, pathway_name)
          kegg_idx <- kegg_idx + 1
        }
      }

      # Trim any unused rows
      if (kegg_idx <= kegg_count) {
        kegg_refs <- kegg_refs[1:(kegg_idx-1), , drop = FALSE]
      }

      result$kegg_refs <- kegg_refs
    }
  }

  if (debug) {
    message("Extracted ", nrow(result$go_terms), " GO terms and ",
            nrow(result$kegg_refs), " KEGG references")
  }

  return(result)
}

#' Create an empty annotation structure
#'
#' @param accession The accession number
#' @return A list with empty annotation structures
#' @export
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
#' This function stores UniProt annotation in the database or adds GO terms and KEGG references
#' to existing annotations if they're missing.
#'
#' @param con A database connection object.
#' @param blast_result_id The ID of the BLAST result to annotate.
#' @param uniprot_info A list containing UniProt information as returned by extract_uniprot_info.
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE.
#' @param update_existing Logical. If TRUE, add GO terms and KEGG references to existing annotations. Default is TRUE.
#'
#' @return The ID of the annotation.
#'
#' @export
store_annotation <- function(con, blast_result_id, uniprot_info, verbose = TRUE, update_existing = TRUE) {
  if (is.null(uniprot_info)) {
    if (verbose) message("Cannot store NULL annotation")
    return(NULL)
  }

  # Add detailed debugging
  if (verbose) {
    message("Storing annotation for blast_result_id: ", blast_result_id)
    message("UniProt info: accession=", uniprot_info$accession,
            ", entry_name=", ifelse(is.null(uniprot_info$entry_name) || is.na(uniprot_info$entry_name), "NA", uniprot_info$entry_name),
            ", gene_names=", ifelse(is.null(uniprot_info$gene_names) || is.na(uniprot_info$gene_names), "", uniprot_info$gene_names))

    if (!is.null(uniprot_info$go_terms) && is.data.frame(uniprot_info$go_terms)) {
      message("GO terms: ", nrow(uniprot_info$go_terms), " entries")
      if (nrow(uniprot_info$go_terms) > 0) {
        message("First GO term: ", uniprot_info$go_terms$go_id[1], " - ", uniprot_info$go_terms$go_term[1])
      }
    } else {
      message("GO terms: NULL or not a data frame")
    }

    if (!is.null(uniprot_info$kegg_refs) && is.data.frame(uniprot_info$kegg_refs)) {
      message("KEGG refs: ", nrow(uniprot_info$kegg_refs), " entries")
      if (nrow(uniprot_info$kegg_refs) > 0) {
        message("First KEGG ref: ", uniprot_info$kegg_refs$kegg_id[1])
      }
    } else {
      message("KEGG refs: NULL or not a data frame")
    }
  }

  # Check if annotation already exists
  existing <- DBI::dbGetQuery(
    con,
    "SELECT annotation_id FROM annotations
     WHERE blast_result_id = ? AND uniprot_accession = ?",
    params = list(blast_result_id, uniprot_info$accession)
  )

  annotation_id <- NULL

  if (nrow(existing) > 0) {
    annotation_id <- existing$annotation_id[1]
    if (verbose) message("Annotation already exists with ID ", annotation_id)

    # If update_existing is TRUE, check if GO terms and KEGG refs exist
    if (update_existing) {
      # Check if GO terms exist for this annotation
      go_count <- DBI::dbGetQuery(
        con,
        "SELECT COUNT(*) AS count FROM go_terms WHERE annotation_id = ?",
        params = list(annotation_id)
      )$count

      # Check if KEGG refs exist for this annotation
      kegg_count <- DBI::dbGetQuery(
        con,
        "SELECT COUNT(*) AS count FROM kegg_references WHERE annotation_id = ?",
        params = list(annotation_id)
      )$count

      if (verbose) {
        message("Existing annotation has ", go_count, " GO terms and ", kegg_count, " KEGG references")
      }

      # If no GO terms and no KEGG refs, add them
      if (go_count == 0 && kegg_count == 0) {
        if (verbose) message("Adding missing GO terms and KEGG references to existing annotation")
      } else {
        if (verbose) message("Skipping - annotation already has GO terms or KEGG references")
        return(annotation_id)
      }
    } else {
      # If not updating existing annotations, just return the ID
      return(annotation_id)
    }
  } else {
    # Create new annotation
    tryCatch({
      # Add annotation
      current_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

      # Prepare parameters with proper NULL/NA handling
      entry_name <- uniprot_info$entry_name
      if (is.null(entry_name) || is.na(entry_name)) entry_name <- NULL

      gene_names <- uniprot_info$gene_names
      if (is.null(gene_names) || is.na(gene_names)) gene_names <- ""

      DBI::dbExecute(
        con,
        "INSERT INTO annotations (blast_result_id, uniprot_accession, entry_name, gene_names, retrieval_date)
         VALUES (?, ?, ?, ?, ?)",
        params = list(
          blast_result_id,
          uniprot_info$accession,
          entry_name,
          gene_names,
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

      if (verbose) message("Created new annotation with ID ", annotation_id)

    }, error = function(e) {
      warning("Error creating annotation: ", e$message)
      return(NULL)
    })
  }

  # Add GO terms
  go_count <- 0
  if (!is.null(annotation_id) && !is.null(uniprot_info$go_terms) &&
      is.data.frame(uniprot_info$go_terms) && nrow(uniprot_info$go_terms) > 0) {
    if (verbose) message("Adding ", nrow(uniprot_info$go_terms), " GO terms")

    tryCatch({
      for (i in 1:nrow(uniprot_info$go_terms)) {
        # Get values with proper NULL handling
        go_id <- as.character(uniprot_info$go_terms$go_id[i])
        if (is.na(go_id)) next

        go_term <- as.character(uniprot_info$go_terms$go_term[i])
        if (is.na(go_term)) go_term <- NULL

        go_category <- as.character(uniprot_info$go_terms$go_category[i])
        if (is.na(go_category)) go_category <- NULL

        go_evidence <- as.character(uniprot_info$go_terms$go_evidence[i])
        if (is.na(go_evidence)) go_evidence <- NULL

        if (verbose && i <= 3) message("  Inserting GO term: ", go_id, " - ", ifelse(is.null(go_term), "NULL", go_term))

        # Execute the insert for this GO term
        result <- DBI::dbExecute(
          con,
          "INSERT INTO go_terms (annotation_id, go_id, go_term, go_category, go_evidence)
           VALUES (?, ?, ?, ?, ?)",
          params = list(
            annotation_id,
            go_id,
            go_term,
            go_category,
            go_evidence
          )
        )

        go_count <- go_count + result
      }

      if (verbose) message("Successfully inserted ", go_count, " GO terms")

    }, error = function(e) {
      warning("Error inserting GO terms: ", e$message)
    })
  } else if (verbose) {
    message("No GO terms to add")
  }

  # Add KEGG references
  kegg_count <- 0
  if (!is.null(annotation_id) && !is.null(uniprot_info$kegg_refs) &&
      is.data.frame(uniprot_info$kegg_refs) && nrow(uniprot_info$kegg_refs) > 0) {
    if (verbose) message("Adding ", nrow(uniprot_info$kegg_refs), " KEGG references")

    tryCatch({
      for (i in 1:nrow(uniprot_info$kegg_refs)) {
        # Get values with proper NULL handling
        kegg_id <- as.character(uniprot_info$kegg_refs$kegg_id[i])
        if (is.na(kegg_id)) next

        pathway_name <- as.character(uniprot_info$kegg_refs$pathway_name[i])
        if (is.na(pathway_name)) pathway_name <- NULL

        if (verbose && i <= 3) message("  Inserting KEGG reference: ", kegg_id)

        # Execute the insert for this KEGG reference
        result <- DBI::dbExecute(
          con,
          "INSERT INTO kegg_references (annotation_id, kegg_id, pathway_name)
           VALUES (?, ?, ?)",
          params = list(
            annotation_id,
            kegg_id,
            pathway_name
          )
        )

        kegg_count <- kegg_count + result
      }

      if (verbose) message("Successfully inserted ", kegg_count, " KEGG references")

    }, error = function(e) {
      warning("Error inserting KEGG references: ", e$message)
    })
  } else if (verbose) {
    message("No KEGG references to add")
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
#' @param delay Delay between operations in seconds. Default is 1.
#' @param offline_mode If TRUE, skips UniProt API and uses basic annotations. Default is FALSE.
#' @param use_cache If TRUE, uses cached API responses if available. Default is TRUE.
#' @param store_cache If TRUE, stores API responses in the database. Default is TRUE.
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE.
#' @param debug_accessions Vector of accessions to debug. Default is NULL.
#'
#' @return A list containing annotation statistics.
#'
#' @export
annotate_blast_results <- function(con, blast_param_id, max_hits = 5, e_value_threshold = 1e-10,
                                   batch_size = 500, delay = 1, offline_mode = FALSE,
                                   use_cache = TRUE, store_cache = TRUE, verbose = TRUE,
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

    # Count existing entries
    if (verbose) {
      cache_count <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS count FROM uniprot_cache")$count
      message("Current UniProt cache contains ", cache_count, " entries")
    }
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

  # Find cached JSON files
  json_dir <- "uniprot_debug"
  cached_accessions <- character(0)
  if (dir.exists(json_dir)) {
    json_files <- list.files(json_dir, pattern = "\\.json$", full.names = FALSE)
    json_files <- grep("_raw\\.json$", json_files, value = TRUE, invert = TRUE)

    if (verbose) message("Found ", length(json_files), " cached JSON files")

    # Extract accessions from filenames
    cached_accessions <- gsub("^uniprot_(.+)\\.json$", "\\1", json_files)
    if (verbose && length(cached_accessions) > 0) {
      message("Cached accessions: ", paste(cached_accessions, collapse = ", "))
    }
  }

  # Debug mode announcement
  if (!is.null(debug_accessions) && length(debug_accessions) > 0) {
    if (verbose) message("Debug mode enabled for accessions: ", paste(debug_accessions, collapse = ", "))
    if (verbose && !offline_mode) message("API queries will still be performed for accessions not in debug mode or cache")
  }

  # Test UniProt API connection only if not in offline mode and not exclusively using cache
  need_api <- !offline_mode &&
    (length(setdiff(hit_accessions, cached_accessions)) > 0) &&
    (!use_cache || store_cache)

  if (need_api) {
    if (verbose) message("Testing UniProt API connection...")
    check_uniprot_connection(verbose = verbose)
  } else {
    if (offline_mode) {
      if (verbose) message("Offline mode enabled - skipping API connection test")
    } else if (use_cache && !store_cache && length(cached_accessions) > 0) {
      if (verbose) message("Using cache only - skipping API connection test")
    }
  }

  # Process all accessions
  uniprot_data <- list()
  cache_hits <- 0
  api_calls <- 0
  cache_storage_attempts <- 0
  cache_storage_success <- 0

  for (i in seq_along(hit_accessions)) {
    acc <- hit_accessions[i]
    enable_debug <- !is.null(debug_accessions) && acc %in% debug_accessions

    # First check if we can get this from the cache directory
    if (acc %in% cached_accessions) {
      if (verbose && enable_debug) message("Reading cached JSON file for ", acc)
      info <- read_uniprot_json(acc, debug = enable_debug)
      uniprot_data[[acc]] <- info
      cache_hits <- cache_hits + 1
    }
    # Then check if we can get this from the database cache
    else if (use_cache) {
      cache_data <- get_cached_uniprot_data(con, acc)
      if (!is.null(cache_data)) {
        if (verbose && enable_debug) message("Retrieved ", acc, " from database cache")
        uniprot_data[[acc]] <- extract_uniprot_info(cache_data, debug = enable_debug)
        cache_hits <- cache_hits + 1
      }
      # Only query API if not in offline mode and not found in any cache
      else if (!offline_mode) {
        if (verbose && enable_debug) message("Querying API for ", acc)
        # Add delay if needed
        if (api_calls > 0 && delay > 0) {
          Sys.sleep(delay)
        }

        # Use the new direct API function
        info <- fetch_uniprot_data(acc, debug = enable_debug)
        api_calls <- api_calls + 1

        # Store the response directly in memory
        uniprot_data[[acc]] <- info

        # Store in database cache if requested
        if (store_cache) {
          cache_storage_attempts <- cache_storage_attempts + 1
          # Create a simplified response to store
          json_data <- jsonlite::toJSON(info, auto_unbox = TRUE)

          tryCatch({
            current_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
            insert_query <- "INSERT INTO uniprot_cache (accession, response_json, retrieval_date) VALUES (?, ?, ?)"
            DBI::dbExecute(con, insert_query, params = list(acc, json_data, current_time))
            cache_storage_success <- cache_storage_success + 1

            if (enable_debug) message("Successfully stored ", acc, " in cache")
          }, error = function(e) {
            if (verbose && enable_debug) {
              message("Error storing ", acc, " in cache: ", e$message)
            }
          })
        }
      }
      else {
        # In offline mode, create empty info
        uniprot_data[[acc]] <- create_empty_annotation(acc)
      }
    }
    # Direct API query if not using cache
    else if (!offline_mode) {
      if (verbose && enable_debug) message("Querying API for ", acc)
      # Add delay if needed
      if (api_calls > 0 && delay > 0) {
        Sys.sleep(delay)
      }

      # Use the new direct API function
      info <- fetch_uniprot_data(acc, debug = enable_debug)
      api_calls <- api_calls + 1

      # Store the response directly in memory
      uniprot_data[[acc]] <- info

      # Store in database cache if requested
      if (store_cache) {
        cache_storage_attempts <- cache_storage_attempts + 1
        # Create a simplified response to store
        json_data <- jsonlite::toJSON(info, auto_unbox = TRUE)

        tryCatch({
          current_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
          insert_query <- "INSERT INTO uniprot_cache (accession, response_json, retrieval_date) VALUES (?, ?, ?)"
          DBI::dbExecute(con, insert_query, params = list(acc, json_data, current_time))
          cache_storage_success <- cache_storage_success + 1

          if (enable_debug) message("Successfully stored ", acc, " in cache")
        }, error = function(e) {
          if (verbose && enable_debug) {
            message("Error storing ", acc, " in cache: ", e$message)
          }
        })
      }
    }
    else {
      # In offline mode, create empty info
      uniprot_data[[acc]] <- create_empty_annotation(acc)
    }

    # Show progress
    if (verbose && (i %% 10 == 0 || i == length(hit_accessions))) {
      message("Processed ", i, " of ", length(hit_accessions), " accessions (",
              cache_hits, " from cache, ", api_calls, " API calls)")
    }

    # Add delay between operations if needed
    if (i < length(hit_accessions) && delay > 0 && api_calls > 0) {
      Sys.sleep(delay)
    }
  }

  if (verbose) {
    message("UniProt data retrieval complete: ", cache_hits, " from cache, ",
            api_calls, " from API")
    if (store_cache) {
      message("Cache storage stats: ", cache_storage_success, " successes out of ",
              cache_storage_attempts, " attempts")
    }
  }

  # No need for separate extraction since we're using direct annotation objects
  if (verbose) message("Successfully extracted information for ", length(uniprot_data), " of ", length(hit_accessions), " accessions.")

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

        # Store annotation using the data we already have
        annotation_id <- store_annotation(
          con,
          blast_result_id,
          uniprot_data[[hit_accession]],
          verbose = (!is.null(debug_accessions) && hit_accession %in% debug_accessions),
          update_existing = TRUE
        )

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

  # Count stored GO terms and KEGG references
  if (verbose) message("Counting stored GO terms and KEGG references...")

  # Get a total count from the database for the specific blast parameter ID
  blast_params_clause <- ""
  if (!is.null(blast_param_id)) {
    blast_params_clause <- paste0(" AND br.blast_param_id = ", blast_param_id)
  }

  go_count <- DBI::dbGetQuery(
    con,
    paste0("SELECT COUNT(*) AS count FROM go_terms g
          JOIN annotations a ON g.annotation_id = a.annotation_id
          JOIN blast_results br ON a.blast_result_id = br.blast_result_id
          WHERE 1=1", blast_params_clause)
  )$count

  kegg_count <- DBI::dbGetQuery(
    con,
    paste0("SELECT COUNT(*) AS count FROM kegg_references k
          JOIN annotations a ON k.annotation_id = a.annotation_id
          JOIN blast_results br ON a.blast_result_id = br.blast_result_id
          WHERE 1=1", blast_params_clause)
  )$count

  if (verbose) {
    message("Stored ", go_count, " GO terms and ", kegg_count, " KEGG references.")

    if (use_cache || store_cache) {
      cache_stats <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS count FROM uniprot_cache")$count
      message("UniProt cache now contains ", cache_stats, " entries")

      # Show new cache entries
      if (cache_stats > 0 && cache_storage_success > 0) {
        message("Successfully created ", cache_storage_success, " new cache entries")
      }
    }

    if (offline_mode) {
      message("\nNOTE: Annotation was performed in offline mode.")
      message("You can re-run the annotation later with UniProt API access to get GO terms and KEGG pathways.")
    }
  }

  # Return summary
  list(
    blast_param_id = blast_param_id,
    unique_accessions = length(hit_accessions),
    successful_extractions = length(uniprot_data),
    annotated_results = successful_annotations,
    go_terms = go_count,
    kegg_refs = kegg_count,
    offline_mode = offline_mode,
    cache_used = use_cache,
    cache_stored = store_cache,
    cache_stored_count = cache_storage_success
  )
}

#' Query the UniProt API for a protein accession
#'
#' @param accession The UniProt accession number
#' @param fields The fields to retrieve from UniProt
#' @param use_rest_api If TRUE, use REST API, otherwise use legacy API
#' @param debug If TRUE, save the API response to a file for debugging and print additional information
#' @param max_retries Maximum number of retry attempts for failed requests
#' @param retry_delay Delay in seconds between retry attempts
#'
#' @return A list containing API response with consistent structure
#' @export
query_uniprot_api <- function(accession,
                              fields = "accession,id,gene_names,go,xref_kegg,organism_name,protein_name",
                              use_rest_api = TRUE,
                              debug = FALSE,
                              max_retries = 3,
                              retry_delay = 2) {

  # Create result structure
  result <- list(
    accession = accession,
    url = NA_character_,
    status_code = NA_integer_,
    content = NA_character_,
    data = NULL,
    error = NA_character_
  )

  # Use API client headers
  headers <- c(
    "User-Agent" = "funseqR/R API client",
    "Accept" = "application/json"
  )

  # Construct URL
  if (use_rest_api) {
    base_url <- "https://rest.uniprot.org/uniprotkb/search"
    query <- paste0("accession:", accession)
    url <- paste0(base_url, "?query=", URLencode(query),
                  "&fields=", URLencode(fields), "&format=json")
  } else {
    base_url <- "https://www.uniprot.org/uniprot"
    query <- paste0("accession:", accession)
    url <- paste0(base_url, "?query=", URLencode(query),
                  "&columns=", URLencode(fields), "&format=json")
  }

  result$url <- url

  # Make HTTP request with retries
  response <- NULL
  success <- FALSE

  for (try_count in 0:max_retries) {
    if (try_count > 0) {
      if (debug) message("Retry attempt ", try_count, " of ", max_retries)
      Sys.sleep(retry_delay)
    }

    tryCatch({
      if (debug) message("Querying URL: ", url)
      response <- httr::GET(url, httr::add_headers(.headers = headers))

      # Get status code
      status <- httr::status_code(response)
      result$status_code <- status

      if (status == 200) {
        if (debug) message("Request successful, status code: 200")
        success <- TRUE
        break  # Exit the retry loop on success
      } else {
        if (debug) message("Request failed with status code: ", status)
      }
    }, error = function(e) {
      if (debug) message("HTTP request error: ", e$message)
      result$error <- paste("HTTP request error:", e$message)
    })

    # Break if last attempt
    if (try_count == max_retries) {
      if (debug) message("All retry attempts failed")
      break
    }
  }

  # If the request failed completely, try fallback or return error
  if (!success) {
    if (use_rest_api) {
      if (debug) message("REST API failed, trying legacy API...")
      return(query_uniprot_api(accession, fields, use_rest_api = FALSE,
                               debug = debug, max_retries = max_retries,
                               retry_delay = retry_delay))
    } else {
      warning(paste("Failed to retrieve data for", accession))
      return(result)
    }
  }

  # Extract and process content
  tryCatch({
    # Get raw content as text
    raw_content <- httr::content(response, "text", encoding = "UTF-8")

    # Important: check if raw_content is not NULL or NA and has content
    if (!is.null(raw_content) && !is.na(raw_content) && nchar(raw_content) > 0) {
      # Save raw content for debugging
      if (debug) {
        debug_dir <- "uniprot_debug"
        if (!dir.exists(debug_dir)) dir.create(debug_dir)
        debug_file <- file.path(debug_dir, paste0("uniprot_", accession, "_raw.json"))
        writeLines(raw_content, debug_file)
        message("Saved raw JSON response to ", debug_file)
      }

      # Update result with content - this was missing in the original function
      result$content <- raw_content

      # Parse JSON
      tryCatch({
        parsed_data <- jsonlite::fromJSON(raw_content, flatten = TRUE)
        result$data <- parsed_data

        if (debug) {
          message("Successfully parsed JSON")
          if (!is.null(parsed_data$results) && length(parsed_data$results) > 0) {
            message("First result has these fields: ",
                    paste(names(parsed_data$results[[1]]), collapse = ", "))
          } else {
            message("No results found in parsed data")
          }
        }
      }, error = function(e) {
        if (debug) message("Error parsing JSON: ", e$message)
        result$error <- paste("Error parsing JSON:", e$message)
      })
    } else {
      if (debug) message("Empty or NULL response content")
      result$error <- "Empty response content"
    }

    return(result)

  }, error = function(e) {
    if (debug) {
      message("Error processing response: ", e$message)
    }

    result$error <- paste("Error processing response:", e$message)
    return(result)
  })

  return(result)  # Final return in case all else fails
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


#' Read UniProt JSON data from a file and extract annotations
#'
#' This function reads a UniProt JSON file and extracts annotation information
#'
#' @param accession The UniProt accession number
#' @param json_file The path to the JSON file (if NULL, will look in uniprot_debug directory)
#' @param debug If TRUE, print debugging information
#'
#' @return A list containing extracted annotation information
#' @export
read_uniprot_json <- function(accession, json_file = NULL, debug = FALSE) {
  # Determine the file path
  if (is.null(json_file)) {
    debug_dir <- "uniprot_debug"
    json_file <- file.path(debug_dir, paste0("uniprot_", accession, ".json"))
  }

  if (!file.exists(json_file)) {
    if (debug) message("JSON file not found: ", json_file)
    return(create_empty_annotation(accession))
  }

  # Read the JSON file
  if (debug) message("Reading JSON file: ", json_file)

  tryCatch({
    # Read as a single string
    json_text <- readChar(json_file, file.info(json_file)$size)
    if (debug) message("Read ", nchar(json_text), " characters")

    # Parse with explicit parameters to prevent simplification
    parsed <- jsonlite::fromJSON(json_text, simplifyVector = FALSE)

    # Create result structure
    result <- list(
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

    # Check for valid structure
    if (is.list(parsed) && "results" %in% names(parsed) &&
        length(parsed$results) > 0 && is.list(parsed$results[[1]])) {

      # Get first result
      entry <- parsed$results[[1]]

      # Basic info
      if ("primaryAccession" %in% names(entry))
        result$accession <- entry$primaryAccession

      if ("uniProtkbId" %in% names(entry))
        result$entry_name <- entry$uniProtkbId

      # Gene names
      if ("genes" %in% names(entry) && length(entry$genes) > 0) {
        gene_names <- c()

        for (i in seq_along(entry$genes)) {
          gene <- entry$genes[[i]]
          if (is.list(gene) && "geneName" %in% names(gene) &&
              is.list(gene$geneName) && "value" %in% names(gene$geneName)) {
            gene_names <- c(gene_names, gene$geneName$value)
          }
        }

        if (length(gene_names) > 0) {
          result$gene_names <- paste(gene_names, collapse = ";")
        }
      }

      # Handle cross-references
      if ("uniProtKBCrossReferences" %in% names(entry) &&
          length(entry$uniProtKBCrossReferences) > 0) {

        # Process GO terms
        go_terms <- data.frame(
          go_id = character(0),
          go_term = character(0),
          go_category = character(0),
          go_evidence = character(0),
          stringsAsFactors = FALSE
        )

        # Process KEGG references
        kegg_refs <- data.frame(
          kegg_id = character(0),
          pathway_name = character(0),
          stringsAsFactors = FALSE
        )

        for (i in seq_along(entry$uniProtKBCrossReferences)) {
          ref <- entry$uniProtKBCrossReferences[[i]]

          # Process GO terms
          if (is.list(ref) && "database" %in% names(ref) &&
              ref$database == "GO" && "id" %in% names(ref)) {

            go_id <- ref$id
            go_term <- NA_character_
            go_evidence <- NA_character_

            if ("properties" %in% names(ref) && length(ref$properties) > 0) {
              for (j in seq_along(ref$properties)) {
                prop <- ref$properties[[j]]
                if (is.list(prop) && "key" %in% names(prop) && "value" %in% names(prop)) {
                  if (prop$key == "GoTerm") go_term <- prop$value
                  if (prop$key == "GoEvidenceType") go_evidence <- prop$value
                }
              }
            }

            if (!is.na(go_term)) {
              new_row <- data.frame(
                go_id = go_id,
                go_term = go_term,
                go_category = substr(go_term, 1, 1),
                go_evidence = go_evidence,
                stringsAsFactors = FALSE
              )
              go_terms <- rbind(go_terms, new_row)
            }
          }

          # Process KEGG references
          if (is.list(ref) && "database" %in% names(ref) &&
              ref$database == "KEGG" && "id" %in% names(ref)) {

            kegg_id <- ref$id
            pathway_name <- NA_character_

            if ("properties" %in% names(ref) && length(ref$properties) > 0) {
              for (j in seq_along(ref$properties)) {
                prop <- ref$properties[[j]]
                if (is.list(prop) && "key" %in% names(prop) && "value" %in% names(prop)) {
                  if (prop$key == "Description") pathway_name <- prop$value
                }
              }
            }

            new_row <- data.frame(
              kegg_id = kegg_id,
              pathway_name = pathway_name,
              stringsAsFactors = FALSE
            )
            kegg_refs <- rbind(kegg_refs, new_row)
          }
        }

        result$go_terms <- go_terms
        result$kegg_refs <- kegg_refs
      }
    }

    if (debug) {
      message("Extracted ", nrow(result$go_terms), " GO terms and ",
              nrow(result$kegg_refs), " KEGG references")
    }

    return(result)

  }, error = function(e) {
    if (debug) message("Error parsing JSON file: ", e$message)
    return(create_empty_annotation(accession))
  })
}

#' Fetch and Process UniProt Data for a Protein Accession
#'
#' This function queries the UniProt REST API for protein information,
#' extracts relevant data, and saves the JSON response to a file for caching.
#' It directly returns a structured annotation object ready for database storage.
#'
#' @param accession The UniProt accession number to query
#' @param debug If TRUE, print detailed debugging information
#' @param timeout Timeout in seconds for the HTTP request. Default is 10.
#' @param retry_max Maximum number of retry attempts. Default is 3.
#' @param retry_delay Delay in seconds between retries. Default is 2.
#'
#' @return A list containing extracted protein information
#' @export
fetch_uniprot_data <- function(accession, debug = FALSE, timeout = 10,
                               retry_max = 3, retry_delay = 2) {
  # Construct URL
  base_url <- "https://rest.uniprot.org/uniprotkb/search"
  fields <- "accession,id,gene_names,go,xref_kegg,organism_name,protein_name"
  query <- paste0("accession:", accession)
  url <- paste0(base_url, "?query=", URLencode(query),
                "&fields=", URLencode(fields), "&format=json")

  if (debug) message("Fetching data from: ", url)

  # Create directory for JSON files if it doesn't exist
  debug_dir <- "uniprot_debug"
  if (!dir.exists(debug_dir)) dir.create(debug_dir)
  json_file <- file.path(debug_dir, paste0("uniprot_", accession, ".json"))

  # Try direct download with download.file or curl
  tryCatch({
    # First try a simpler approach with download.file
    if (debug) message("Using download.file to fetch data")
    temp_file <- tempfile()
    download.file(url, temp_file, quiet = !debug, mode = "wb",
                  headers = c("Accept" = "application/json",
                              "User-Agent" = "Mozilla/5.0 (R; funseqR)"))

    # Check if the file has content
    if (file.size(temp_file) > 10) {
      # Read the file
      content_text <- readChar(temp_file, file.size(temp_file))

      if (debug) {
        message("Downloaded content length: ", nchar(content_text))
        message("First 100 chars: ", substr(content_text, 1, 100))
      }

      # Try to parse to validate
      tryCatch({
        parsed <- jsonlite::fromJSON(content_text, flatten = TRUE)

        # If we get here, the JSON is valid
        # Write to the cache file
        writeLines(content_text, json_file)

        if (debug) message("Valid JSON saved to: ", json_file)

        # Extract and return the annotations
        return(process_uniprot_json(parsed, accession, debug))
      }, error = function(e) {
        if (debug) message("Error parsing downloaded JSON: ", e$message)
      })
    } else {
      if (debug) message("Downloaded file is empty or too small")
    }

    # If we get here, the download.file approach failed, try with httr
    if (debug) message("Trying alternative approach with httr")

    # Use httr with explicit timeouts
    response <- httr::GET(
      url,
      httr::add_headers(
        "Accept" = "application/json",
        "User-Agent" = "Mozilla/5.0 (R; funseqR)"
      ),
      httr::timeout(timeout)
    )

    status <- httr::status_code(response)

    if (debug) message("HTTP status code: ", status)

    if (status == 200) {
      # Get content as text with explicit encoding
      content_text <- httr::content(response, "text", encoding = "UTF-8")

      if (debug) {
        message("Content length: ", nchar(content_text))
        message("First 100 chars: ", substr(content_text, 1, 100))
      }

      # Make sure we got valid content
      if (!is.null(content_text) && nchar(content_text) > 10) {
        # Try to parse to validate
        tryCatch({
          parsed <- jsonlite::fromJSON(content_text, flatten = TRUE)

          # Write to the cache file
          writeLines(content_text, json_file)

          if (debug) message("Valid JSON saved to: ", json_file)

          # Extract and return the annotations
          return(process_uniprot_json(parsed, accession, debug))
        }, error = function(e) {
          if (debug) message("Error parsing JSON: ", e$message)
          # Save the raw content even if parsing failed
          writeLines(content_text, json_file)
          if (debug) message("Raw content saved to: ", json_file)
        })
      } else {
        if (debug) message("Empty or invalid content received")
        writeLines("Empty or invalid content", json_file)
      }
    } else {
      if (debug) message("HTTP request failed with status: ", status)
      writeLines(paste("HTTP error:", status), json_file)
    }
  }, error = function(e) {
    if (debug) message("Error fetching data: ", e$message)
    writeLines(paste("Error:", e$message), json_file)
  })

  # If we reach here, the fetch failed - create empty annotation
  if (debug) message("Returning empty annotation for ", accession)
  return(create_empty_annotation(accession))
}

#' Process UniProt JSON data into annotations
#'
#' Helper function to extract annotation information from parsed UniProt JSON
#'
#' @param parsed Parsed JSON data from UniProt API
#' @param accession The accession number
#' @param debug Enable debug output
#'
#' @return A structured annotation object
process_uniprot_json <- function(parsed, accession, debug = FALSE) {
  # Create empty annotation structure
  info <- create_empty_annotation(accession)

  # Process only if we have results
  if (!is.null(parsed$results) && length(parsed$results) > 0) {
    entry <- parsed$results[[1]]

    # Basic info
    if (!is.null(entry$primaryAccession))
      info$accession <- entry$primaryAccession

    if (!is.null(entry$uniProtkbId))
      info$entry_name <- entry$uniProtkbId

    # Gene names
    if (!is.null(entry$genes) && length(entry$genes) > 0) {
      gene_names <- c()

      for (gene in entry$genes) {
        if (!is.null(gene$geneName) && !is.null(gene$geneName$value)) {
          gene_names <- c(gene_names, gene$geneName$value)
        }
      }

      if (length(gene_names) > 0) {
        info$gene_names <- paste(gene_names, collapse = ";")
      }
    }

    # GO terms and KEGG references
    if (!is.null(entry$uniProtKBCrossReferences) && length(entry$uniProtKBCrossReferences) > 0) {
      go_terms <- data.frame(
        go_id = character(0),
        go_term = character(0),
        go_category = character(0),
        go_evidence = character(0),
        stringsAsFactors = FALSE
      )

      kegg_refs <- data.frame(
        kegg_id = character(0),
        pathway_name = character(0),
        stringsAsFactors = FALSE
      )

      for (ref in entry$uniProtKBCrossReferences) {
        # GO terms
        if (!is.null(ref$database) && ref$database == "GO" && !is.null(ref$id)) {
          go_id <- ref$id
          go_term <- NA_character_
          go_evidence <- NA_character_

          if (!is.null(ref$properties)) {
            for (prop in ref$properties) {
              if (!is.null(prop$key) && !is.null(prop$value)) {
                if (prop$key == "GoTerm") go_term <- prop$value
                if (prop$key == "GoEvidenceType") go_evidence <- prop$value
              }
            }
          }

          if (!is.na(go_term)) {
            new_row <- data.frame(
              go_id = go_id,
              go_term = go_term,
              go_category = substr(go_term, 1, 1),
              go_evidence = go_evidence,
              stringsAsFactors = FALSE
            )
            go_terms <- rbind(go_terms, new_row)
          }
        }

        # KEGG references
        if (!is.null(ref$database) && ref$database == "KEGG" && !is.null(ref$id)) {
          kegg_id <- ref$id
          pathway_name <- NA_character_

          if (!is.null(ref$properties)) {
            for (prop in ref$properties) {
              if (!is.null(prop$key) && !is.null(prop$value)) {
                if (prop$key == "Description") pathway_name <- prop$value
              }
            }
          }

          new_row <- data.frame(
            kegg_id = kegg_id,
            pathway_name = pathway_name,
            stringsAsFactors = FALSE
          )
          kegg_refs <- rbind(kegg_refs, new_row)
        }
      }

      info$go_terms <- go_terms
      info$kegg_refs <- kegg_refs
    }
  }

  if (debug) {
    message("Extracted ", nrow(info$go_terms), " GO terms and ",
            nrow(info$kegg_refs), " KEGG references")
  }

  return(info)
}
