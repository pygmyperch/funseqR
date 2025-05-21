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
#' @param debug If TRUE, print debugging information
#' @return The cached UniProt data or NULL if not found
#' @export
get_cached_uniprot_data <- function(con, accession, debug = FALSE) {
  if (is.null(con) || !DBI::dbIsValid(con)) {
    if (debug) message("Invalid database connection")
    return(NULL)
  }

  # Check if the table exists
  tables <- DBI::dbListTables(con)
  if (!"uniprot_cache" %in% tables) {
    if (debug) message("UniProt cache table doesn't exist")
    return(NULL)
  }

  # Check if data exists in cache
  tryCatch({
    cache_query <- "SELECT response_json, retrieval_date FROM uniprot_cache WHERE accession = ?"
    cache_result <- DBI::dbGetQuery(con, cache_query, params = list(accession))

    if (nrow(cache_result) == 0 || is.null(cache_result$response_json)) {
      if (debug) message("No cache entry found for ", accession)
      return(NULL)
    }

    # Get the JSON text
    json_text <- cache_result$response_json[1]
    retrieval_date <- cache_result$retrieval_date[1]

    if (is.na(json_text) || nchar(json_text) == 0) {
      if (debug) message("Cache entry for ", accession, " is empty")
      return(NULL)
    }

    if (debug) {
      message("Found cache entry for ", accession, " retrieved on ", retrieval_date)
      message("JSON length: ", nchar(json_text), " chars")
      if (nchar(json_text) > 0) {
        message("First 50 chars: ", substr(json_text, 1, 50), "...")
      }
    }

    # Check for placeholder entries
    if (grepl('"empty"\\s*:\\s*true', json_text)) {
      if (debug) message("Cache entry is a placeholder, ignoring")
      return(NULL)
    }

    # Create a response object - without trying to parse the JSON
    # We'll leave the parsing to extract_uniprot_info to avoid duplicate parsing
    result <- list(
      accession = accession,
      url = NA_character_,
      status_code = 200,
      content = json_text,   # Original JSON text
      data = NULL,           # We won't parse here to avoid issues
      error = NA_character_,
      from_cache = TRUE,
      retrieval_date = retrieval_date
    )

    return(result)
  }, error = function(e) {
    if (debug) message("Error retrieving cache data for ", accession, ": ", e$message)
    return(NULL)
  })
}

#' Store UniProt data in cache with robust error handling
#'
#' @param con A database connection object
#' @param accession The UniProt accession number
#' @param response The UniProt API response
#' @param debug If TRUE, print debugging information
#' @return TRUE if storage successful, FALSE otherwise
#' @export
store_uniprot_data <- function(con, accession, response, debug = FALSE) {
  # Initialize json_data to a fallback value
  json_data <- paste0('{"accession":"', accession, '","empty":true,"timestamp":"',
                      format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '"}')

  if (debug) {
    message("Attempting to store data for accession: ", accession)
  }

  # First check if the connection is valid
  if (is.null(con) || !DBI::dbIsValid(con)) {
    if (debug) message("Database connection is invalid or NULL")
    return(FALSE)
  }

  # Handle the response object carefully with separate conditionals
  if (!is.null(response)) {
    # Determine if it's a list
    if (is.list(response)) {
      # Try content field if it exists
      if (!is.null(names(response)) && "content" %in% names(response)) {
        # Content exists, check if it's usable
        content_val <- response$content
        if (!is.null(content_val) && is.character(content_val) && !is.na(content_val) && nchar(content_val) > 0) {
          json_data <- content_val
          if (debug) message("Using content field as JSON")
        }
      }
      # If content didn't work, try data field
      else if (!is.null(names(response)) && "data" %in% names(response) && !is.null(response$data)) {
        tryCatch({
          json_text <- jsonlite::toJSON(response$data, auto_unbox = TRUE)
          if (!is.null(json_text) && nchar(json_text) > 0) {
            json_data <- json_text
            if (debug) message("Converted data field to JSON")
          }
        }, error = function(e) {
          if (debug) message("Error converting data to JSON: ", e$message)
        })
      }
      # As last resort, try to convert the whole response
      else {
        tryCatch({
          json_text <- jsonlite::toJSON(response, auto_unbox = TRUE)
          if (!is.null(json_text) && nchar(json_text) > 0) {
            json_data <- json_text
            if (debug) message("Converted entire response to JSON")
          }
        }, error = function(e) {
          if (debug) message("Error converting response to JSON: ", e$message)
        })
      }
    }
    else if (is.character(response) && !is.na(response) && nchar(response) > 0) {
      # The response itself is a string, try to use it directly
      json_data <- response
      if (debug) message("Using response as JSON string")
    }
  }

  # Store in database
  success <- FALSE
  tryCatch({
    current_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

    # First check if we need to create the table
    tables <- DBI::dbListTables(con)
    if (!"uniprot_cache" %in% tables) {
      if (debug) message("Creating uniprot_cache table")
      DBI::dbExecute(con, "
        CREATE TABLE uniprot_cache (
          cache_id INTEGER PRIMARY KEY,
          accession TEXT NOT NULL UNIQUE,
          response_json TEXT NOT NULL,
          retrieval_date TEXT NOT NULL
        )
      ")
      DBI::dbExecute(con, "CREATE INDEX idx_uniprot_cache_accession ON uniprot_cache (accession)")
    }

    # Check if entry exists
    existing_query <- "SELECT COUNT(*) AS count FROM uniprot_cache WHERE accession = ?"
    existing <- DBI::dbGetQuery(con, existing_query, params = list(accession))
    has_existing <- FALSE
    if (nrow(existing) > 0 && !is.null(existing$count) && existing$count > 0) {
      has_existing <- TRUE
    }

    if (has_existing) {
      # Update existing record
      update_query <- "UPDATE uniprot_cache SET response_json = ?, retrieval_date = ? WHERE accession = ?"
      result <- DBI::dbExecute(con, update_query, params = list(json_data, current_time, accession))
      if (debug) message("Updated existing record for ", accession)
      success <- TRUE
    } else {
      # Insert new record
      insert_query <- "INSERT INTO uniprot_cache (accession, response_json, retrieval_date) VALUES (?, ?, ?)"
      result <- DBI::dbExecute(con, insert_query, params = list(accession, json_data, current_time))
      if (debug) message("Created new record for ", accession)
      success <- TRUE
    }
  }, error = function(e) {
    if (debug) message("Database error: ", e$message)
  })

  return(success)
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
  # Initialize result with empty data
  result <- list(
    accession = if (!is.null(uniprot_data$accession)) uniprot_data$accession else "Unknown",
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

  if (is.null(uniprot_data) || !is.list(uniprot_data)) {
    if (debug) message("No valid UniProt data provided")
    return(result)
  }

  # Debug the structure of what we received
  if (debug) {
    message("Input data structure:")
    message("Class: ", paste(class(uniprot_data), collapse = ", "))

    if (is.list(uniprot_data)) {
      message("List with fields: ", paste(names(uniprot_data), collapse = ", "))

      if ("content" %in% names(uniprot_data)) {
        if (is.character(uniprot_data$content)) {
          message("Content is character with length ", nchar(uniprot_data$content))
        } else {
          message("Content is not a character but a ", class(uniprot_data$content)[1])
        }
      }

      if ("data" %in% names(uniprot_data)) {
        message("Data is class ", class(uniprot_data$data)[1])
        if (is.list(uniprot_data$data)) {
          message("Data fields: ", paste(names(uniprot_data$data), collapse = ", "))
        }
      }
    }
  }

  # Get the entry data - handle both formats (direct and search)
  entry <- NULL

  # Try different paths to find the entry data

  # Path 1: From content field containing raw JSON - preferred approach
  if (is.null(entry) &&
      !is.null(uniprot_data$content) &&
      is.character(uniprot_data$content) &&
      nchar(uniprot_data$content) > 0) {

    if (debug) message("Trying to extract from content field (raw JSON)")

    tryCatch({
      # Parse the JSON content
      parsed_content <- jsonlite::fromJSON(uniprot_data$content, flatten = FALSE)

      if (debug) {
        message("Successfully parsed JSON from content")
        message("Parsed type: ", class(parsed_content)[1])
        if (is.list(parsed_content)) {
          message("Parsed fields: ", paste(names(parsed_content), collapse = ", "))
        }
      }

      # Check which format we have
      if (is.list(parsed_content)) {
        # Case 1: Search API response with results array
        if (!is.null(parsed_content$results) &&
            is.list(parsed_content$results) &&
            length(parsed_content$results) > 0) {

          if (debug) message("Found results array in content")
          entry <- parsed_content$results[[1]]  # Take first result
        }
        # Case 2: Direct API response with entry data at root
        else if (!is.null(parsed_content$primaryAccession)) {
          if (debug) message("Found direct entry in content")
          entry <- parsed_content
        }
        # Case 3: Unknown format but is a list
        else {
          if (debug) message("Unknown list format in content, using as-is")
          entry <- parsed_content
        }
      }
      # Atomic response - try to use as-is
      else if (length(parsed_content) > 0) {
        if (debug) message("Parsed content is atomic, using as-is")
        # Create a minimal entry with the accession
        entry <- list(primaryAccession = result$accession)
      }
    }, error = function(e) {
      if (debug) message("Error parsing content: ", e$message)
    })
  }

  # Path 2: From data field containing parsed JSON
  if (is.null(entry) && !is.null(uniprot_data$data)) {
    if (debug) message("Trying to extract from data field (parsed object)")

    # Already parsed data
    parsed_data <- uniprot_data$data

    # Handle different types
    if (is.list(parsed_data)) {
      # Case 1: Search API response with results array
      if (!is.null(parsed_data$results) &&
          is.list(parsed_data$results) &&
          length(parsed_data$results) > 0) {

        if (debug) message("Found results array in data")
        entry <- parsed_data$results[[1]]  # Take first result
      }
      # Case 2: Direct API response with entry data at root
      else if (!is.null(parsed_data$primaryAccession)) {
        if (debug) message("Found direct entry in data")
        entry <- parsed_data
      }
      # Case 3: Unknown format but is a list
      else {
        if (debug) message("Unknown list format in data, using as-is")
        entry <- parsed_data
      }
    }
    # Handle atomic data
    else if (length(parsed_data) > 0) {
      if (debug) message("Data is atomic, using as-is")
      # Create a minimal entry with the accession
      entry <- list(primaryAccession = result$accession)
    }
  }

  # Path 3: Check if uniprot_data itself is the entry
  if (is.null(entry) && is.list(uniprot_data) && !is.null(uniprot_data$primaryAccession)) {
    if (debug) message("Found entry data at root level")
    entry <- uniprot_data
  }

  # If we couldn't find entry data, return empty result
  if (is.null(entry)) {
    if (debug) message("No entry data found in response")
    return(result)
  }

  # Make sure entry is a list before trying to access fields
  if (!is.list(entry)) {
    if (debug) message("Entry is not a list, creating minimal entry")
    entry <- list(primaryAccession = result$accession)
  }

  # Extract basic info - with type checking
  if (is.list(entry) && !is.null(entry$primaryAccession)) {
    result$accession <- entry$primaryAccession
    if (debug) message("Found accession: ", result$accession)
  }

  if (is.list(entry) && !is.null(entry$uniProtkbId)) {
    result$entry_name <- entry$uniProtkbId
    if (debug) message("Found entry name: ", result$entry_name)
  }

  # Extract gene names with type checking
  if (is.list(entry) && !is.null(entry$genes) && is.list(entry$genes) && length(entry$genes) > 0) {
    gene_names <- c()

    for (i in seq_along(entry$genes)) {
      gene <- entry$genes[[i]]
      if (is.list(gene) && !is.null(gene$geneName) && is.list(gene$geneName) && !is.null(gene$geneName$value)) {
        gene_names <- c(gene_names, gene$geneName$value)
        if (debug) message("Found gene name: ", gene$geneName$value)
      }
    }

    if (length(gene_names) > 0) {
      result$gene_names <- paste(gene_names, collapse = ";")
      if (debug) message("Combined gene names: ", result$gene_names)
    }
  }

  # Extract GO terms and KEGG references with type checking
  if (is.list(entry) && !is.null(entry$uniProtKBCrossReferences) &&
      is.list(entry$uniProtKBCrossReferences) && length(entry$uniProtKBCrossReferences) > 0) {

    if (debug) {
      message("Processing ", length(entry$uniProtKBCrossReferences), " cross-references")
    }

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

    # Process each cross-reference
    for (i in seq_along(entry$uniProtKBCrossReferences)) {
      ref <- entry$uniProtKBCrossReferences[[i]]

      # Debug the reference type
      if (debug) {
        ref_type <- "unknown"
        if (is.list(ref) && !is.null(ref$database)) {
          ref_type <- ref$database
        }
        message("Reference #", i, " type: ", ref_type)
      }

      # Skip if not a valid reference
      if (!is.list(ref) || is.null(ref$database)) {
        if (debug) message("Skipping invalid reference")
        next
      }

      # Process GO terms
      if (ref$database == "GO" && !is.null(ref$id)) {
        if (debug) message("Processing GO reference: ", ref$id)

        go_id <- ref$id
        go_term <- NA_character_
        go_evidence <- NA_character_

        # Debug the properties structure
        if (debug && !is.null(ref$properties)) {
          message("Properties structure for GO term:")
          print(str(ref$properties))
        }

        # Extract properties
        if (!is.null(ref$properties) && is.list(ref$properties) && length(ref$properties) > 0) {
          for (j in seq_along(ref$properties)) {
            prop <- ref$properties[[j]]

            if (debug) {
              prop_key <- if (!is.null(prop$key)) prop$key else "NULL"
              prop_value <- if (!is.null(prop$value)) prop$value else "NULL"
              message("Property #", j, ": ", prop_key, " = ", prop_value)
            }

            if (is.list(prop) && !is.null(prop$key) && !is.null(prop$value)) {
              if (prop$key == "GoTerm") {
                go_term <- prop$value
                if (debug) message("Found GO term: ", go_term)
              }
              if (prop$key == "GoEvidenceType") {
                go_evidence <- prop$value
                if (debug) message("Found GO evidence: ", go_evidence)
              }
            }
          }
        }

        # Add to GO terms if we have a term
        if (!is.na(go_term)) {
          # Extract category from the term (e.g., "C:nucleus" -> "C")
          go_category <- NA_character_
          if (nchar(go_term) >= 1) {
            go_category <- substr(go_term, 1, 1)
            if (debug) message("Extracted GO category: ", go_category)
          }

          new_row <- data.frame(
            go_id = go_id,
            go_term = go_term,
            go_category = go_category,
            go_evidence = go_evidence,
            stringsAsFactors = FALSE
          )
          go_terms <- rbind(go_terms, new_row)

          if (debug) message("Added GO term: ", go_id, " - ", go_term)
        } else {
          if (debug) message("Skipping GO term (no term value found)")
        }
      }

      # Process KEGG references
      if (ref$database == "KEGG" && !is.null(ref$id)) {
        if (debug) message("Processing KEGG reference: ", ref$id)

        kegg_id <- ref$id
        pathway_name <- NA_character_

        # Debug the properties structure
        if (debug && !is.null(ref$properties)) {
          message("Properties structure for KEGG reference:")
          print(str(ref$properties))
        }

        # Extract properties
        if (!is.null(ref$properties) && is.list(ref$properties) && length(ref$properties) > 0) {
          for (j in seq_along(ref$properties)) {
            prop <- ref$properties[[j]]

            if (debug) {
              prop_key <- if (!is.null(prop$key)) prop$key else "NULL"
              prop_value <- if (!is.null(prop$value)) prop$value else "NULL"
              message("Property #", j, ": ", prop_key, " = ", prop_value)
            }

            if (is.list(prop) && !is.null(prop$key) && !is.null(prop$value)) {
              if (prop$key == "Description" || prop$key == "PathwayName") {
                pathway_name <- prop$value
                if (debug) message("Found pathway name: ", pathway_name)
              }
            }
          }
        }

        # Add KEGG reference
        new_row <- data.frame(
          kegg_id = kegg_id,
          pathway_name = pathway_name,
          stringsAsFactors = FALSE
        )
        kegg_refs <- rbind(kegg_refs, new_row)

        if (debug) message("Added KEGG reference: ", kegg_id)
      }
    }

    # Update result with extracted data
    result$go_terms <- go_terms
    result$kegg_refs <- kegg_refs

    if (debug) {
      message("Extracted ", nrow(go_terms), " GO terms and ",
              nrow(kegg_refs), " KEGG references")
    }
  } else if (debug) {
    message("No cross-references found in entry or structure unexpected")
    # Debug what we're seeing in the entry
    if (is.list(entry)) {
      message("Entry fields: ", paste(names(entry), collapse = ", "))
      if ("uniProtKBCrossReferences" %in% names(entry)) {
        message("uniProtKBCrossReferences type: ", class(entry$uniProtKBCrossReferences)[1])
        message("uniProtKBCrossReferences length: ", length(entry$uniProtKBCrossReferences))
      }
    }
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

  # Check for required fields
  accession <- uniprot_info$accession
  if (is.null(accession) || is.na(accession) || nchar(accession) == 0) {
    if (verbose) message("Missing accession in annotation data")
    return(NULL)
  }

  # Debugging
  if (verbose) {
    message("Storing annotation for blast_result_id: ", blast_result_id)
    message("UniProt info: accession=", accession,
            ", entry_name=", ifelse(is.null(uniprot_info$entry_name) || is.na(uniprot_info$entry_name), "NA", uniprot_info$entry_name),
            ", gene_names=", ifelse(is.null(uniprot_info$gene_names) || is.na(uniprot_info$gene_names), "", uniprot_info$gene_names))

    # Log what we have for GO terms and KEGG refs
    if (!is.null(uniprot_info$go_terms) && is.data.frame(uniprot_info$go_terms)) {
      message("GO terms: ", nrow(uniprot_info$go_terms), " entries")
    } else {
      message("GO terms: NULL or not a data frame")
    }

    if (!is.null(uniprot_info$kegg_refs) && is.data.frame(uniprot_info$kegg_refs)) {
      message("KEGG refs: ", nrow(uniprot_info$kegg_refs), " entries")
    } else {
      message("KEGG refs: NULL or not a data frame")
    }
  }

  # Check if annotation already exists
  existing <- DBI::dbGetQuery(
    con,
    "SELECT annotation_id FROM annotations WHERE blast_result_id = ? AND uniprot_accession = ?",
    params = list(blast_result_id, accession)
  )

  annotation_id <- NULL

  if (nrow(existing) > 0) {
    annotation_id <- existing$annotation_id[1]
    if (verbose) message("Annotation already exists with ID ", annotation_id)

    # Check for GO terms and KEGG refs if updating existing annotations
    if (update_existing) {
      go_count <- DBI::dbGetQuery(
        con,
        "SELECT COUNT(*) AS count FROM go_terms WHERE annotation_id = ?",
        params = list(annotation_id)
      )$count

      kegg_count <- DBI::dbGetQuery(
        con,
        "SELECT COUNT(*) AS count FROM kegg_references WHERE annotation_id = ?",
        params = list(annotation_id)
      )$count

      if (verbose) {
        message("Existing annotation has ", go_count, " GO terms and ", kegg_count, " KEGG references")
      }

      if (go_count == 0 && kegg_count == 0) {
        if (verbose) message("Adding missing GO terms and KEGG references to existing annotation")
      } else {
        if (verbose) message("Skipping - annotation already has GO terms or KEGG references")
        return(annotation_id)
      }
    } else {
      return(annotation_id)
    }
  } else {
    # Create new annotation
    tryCatch({
      # Add annotation
      current_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

      # Handle entry_name properly - if it's NULL, NA, or not length 1, use NULL in the query
      entry_name <- NULL
      if (!is.null(uniprot_info$entry_name) &&
          length(uniprot_info$entry_name) == 1 &&
          !is.na(uniprot_info$entry_name)) {
        entry_name <- uniprot_info$entry_name
      }

      # Handle gene_names
      gene_names <- ""
      if (!is.null(uniprot_info$gene_names) &&
          length(uniprot_info$gene_names) == 1 &&
          !is.na(uniprot_info$gene_names)) {
        gene_names <- uniprot_info$gene_names
      }

      # Use a query that only includes entry_name if it's not NULL
      if (is.null(entry_name)) {
        query <- "INSERT INTO annotations (blast_result_id, uniprot_accession, gene_names, retrieval_date) VALUES (?, ?, ?, ?)"
        params <- list(blast_result_id, accession, gene_names, current_time)
      } else {
        query <- "INSERT INTO annotations (blast_result_id, uniprot_accession, entry_name, gene_names, retrieval_date) VALUES (?, ?, ?, ?, ?)"
        params <- list(blast_result_id, accession, entry_name, gene_names, current_time)
      }

      # Execute the query
      DBI::dbExecute(con, query, params = params)

      # Get the ID of the newly created annotation
      annotation_id <- DBI::dbGetQuery(
        con,
        "SELECT annotation_id FROM annotations WHERE blast_result_id = ? AND uniprot_accession = ? ORDER BY annotation_id DESC LIMIT 1",
        params = list(blast_result_id, accession)
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
        go_id <- uniprot_info$go_terms$go_id[i]
        if (is.na(go_id)) next

        # Check if GO term already exists
        go_exists <- DBI::dbGetQuery(
          con,
          "SELECT COUNT(*) AS count FROM go_terms WHERE annotation_id = ? AND go_id = ?",
          params = list(annotation_id, go_id)
        )$count > 0

        if (go_exists) {
          if (verbose) message("GO term ", go_id, " already exists, skipping")
          next
        }

        # Prepare parameters for insertion
        go_term <- uniprot_info$go_terms$go_term[i]
        go_category <- uniprot_info$go_terms$go_category[i]
        go_evidence <- uniprot_info$go_terms$go_evidence[i]

        # Handle NULL/NA values
        if (is.na(go_term)) go_term <- NULL
        if (is.na(go_category)) go_category <- NULL
        if (is.na(go_evidence)) go_evidence <- NULL

        # Build query based on what's available
        fields <- c("annotation_id", "go_id")
        values <- c("?", "?")
        insert_params <- list(annotation_id, go_id)

        if (!is.null(go_term)) {
          fields <- c(fields, "go_term")
          values <- c(values, "?")
          insert_params <- c(insert_params, list(go_term))
        }

        if (!is.null(go_category)) {
          fields <- c(fields, "go_category")
          values <- c(values, "?")
          insert_params <- c(insert_params, list(go_category))
        }

        if (!is.null(go_evidence)) {
          fields <- c(fields, "go_evidence")
          values <- c(values, "?")
          insert_params <- c(insert_params, list(go_evidence))
        }

        # Create and execute the insertion query
        insert_query <- paste0(
          "INSERT INTO go_terms (", paste(fields, collapse = ", "), ") ",
          "VALUES (", paste(values, collapse = ", "), ")"
        )

        DBI::dbExecute(con, insert_query, params = insert_params)
        go_count <- go_count + 1
      }

      if (verbose) message("Added ", go_count, " GO terms")
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
        kegg_id <- uniprot_info$kegg_refs$kegg_id[i]
        if (is.na(kegg_id)) next

        # Check if KEGG reference already exists
        kegg_exists <- DBI::dbGetQuery(
          con,
          "SELECT COUNT(*) AS count FROM kegg_references WHERE annotation_id = ? AND kegg_id = ?",
          params = list(annotation_id, kegg_id)
        )$count > 0

        if (kegg_exists) {
          if (verbose) message("KEGG reference ", kegg_id, " already exists, skipping")
          next
        }

        # Prepare parameters for insertion
        pathway_name <- uniprot_info$kegg_refs$pathway_name[i]
        if (is.na(pathway_name)) pathway_name <- NULL

        # Build query based on what's available
        if (is.null(pathway_name)) {
          insert_query <- "INSERT INTO kegg_references (annotation_id, kegg_id) VALUES (?, ?)"
          insert_params <- list(annotation_id, kegg_id)
        } else {
          insert_query <- "INSERT INTO kegg_references (annotation_id, kegg_id, pathway_name) VALUES (?, ?, ?)"
          insert_params <- list(annotation_id, kegg_id, pathway_name)
        }

        DBI::dbExecute(con, insert_query, params = insert_params)
        kegg_count <- kegg_count + 1
      }

      if (verbose) message("Added ", kegg_count, " KEGG references")
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
#' @param con Database connection object.
#' @param blast_param_id The ID of the BLAST parameters.
#' @param max_hits Maximum number of hits to process per query. Default is 5.
#' @param e_value_threshold E-value threshold for filtering BLAST hits. Default is 1e-10.
#' @param batch_size Number of annotations to process in each transaction. Default is 500.
#' @param delay Delay between operations in seconds. Default is 1.
#' @param offline_mode If TRUE, skips UniProt API and uses basic annotations. Default is FALSE.
#' @param use_cache If TRUE, uses cached API responses if available. Default is TRUE.
#' @param store_cache If TRUE, stores API responses in the database. Default is TRUE.
#' @param verify_storage If TRUE, verifies annotation storage after processing. Default is TRUE.
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE.
#' @param debug_accessions Vector of accessions to debug. Default is NULL.
#'
#' @return A list containing annotation statistics.
#'
#' @export
annotate_blast_results <- function(con, blast_param_id, max_hits = 5, e_value_threshold = 1e-10,
                                   batch_size = 500, delay = 1, offline_mode = FALSE,
                                   use_cache = TRUE, store_cache = TRUE, verify_storage = TRUE,
                                   verbose = TRUE, debug_accessions = NULL) {
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

  # Test UniProt API connection only if not in offline mode
  if (!offline_mode && !use_cache) {
    if (verbose) message("Testing UniProt API connection...")
    check_uniprot_connection(verbose = verbose)
  } else if (offline_mode) {
    if (verbose) message("Offline mode enabled - skipping API connection test")
  }

  # Process all accessions to get UniProt data
  uniprot_data <- list()
  cache_hits <- 0
  api_calls <- 0
  storage_successes <- 0

  for (i in seq_along(hit_accessions)) {
    acc <- hit_accessions[i]
    enable_debug <- !is.null(debug_accessions) && acc %in% debug_accessions

    # First always check database cache if use_cache is TRUE
    if (use_cache) {
      cache_data <- get_cached_uniprot_data(con, acc)
      if (!is.null(cache_data)) {
        if (verbose && enable_debug) message("Retrieved ", acc, " from database cache")
        uniprot_data[[acc]] <- extract_uniprot_info(cache_data, debug = enable_debug)
        cache_hits <- cache_hits + 1
        next  # Skip to next accession
      }
    }

    # Check JSON files as secondary cache source
    json_file <- file.path("uniprot_debug", paste0("uniprot_", acc, ".json"))
    if (use_cache && file.exists(json_file) && file.size(json_file) > 10) {
      if (verbose && enable_debug) message("Reading ", acc, " from JSON file")

      # Read and parse the file
      tryCatch({
        file_data <- read_uniprot_json(acc, json_file, debug = enable_debug)
        uniprot_data[[acc]] <- file_data

        # Store in database cache if requested
        if (store_cache) {
          # Convert to JSON for storage
          json_text <- readChar(json_file, file.info(json_file)$size)

          # Create a response-like object
          response_obj <- list(
            accession = acc,
            content = json_text,
            status_code = 200
          )

          # Store in database
          storage_success <- store_uniprot_data(con, acc, response_obj, debug = enable_debug)
          if (storage_success) storage_successes <- storage_successes + 1
        }

        cache_hits <- cache_hits + 1
        next  # Skip to next accession
      }, error = function(e) {
        if (verbose && enable_debug) {
          message("Error reading JSON file for ", acc, ": ", e$message)
        }
        # Continue to API fallback
      })
    }

    # If not in cache and not in offline mode, fetch from API
    if (!offline_mode) {
      if (verbose && enable_debug) message("Fetching ", acc, " from UniProt API")

      # Add delay if needed
      if (api_calls > 0 && delay > 0) {
        Sys.sleep(delay)
      }

      # Wrap API call in tryCatch
      api_result <- NULL
      tryCatch({
        api_result <- query_uniprot_api(acc, debug = enable_debug)
        api_calls <- api_calls + 1
      }, error = function(e) {
        if (verbose) message("Error querying API for ", acc, ": ", e$message)
        # Create a minimal result object
        api_result <- list(
          accession = acc,
          error = e$message
        )
      })

      # ALWAYS store API result in database cache if requested
      if (store_cache && !is.null(api_result)) {
        # Add extra debug for storage process
        storage_success <- store_uniprot_data(con, acc, api_result, debug = enable_debug || verbose)
        if (storage_success) storage_successes <- storage_successes + 1
      }

      # Extract information, even if minimal
      info <- NULL
      tryCatch({
        info <- extract_uniprot_info(api_result, debug = enable_debug)
      }, error = function(e) {
        if (verbose) message("Error extracting info for ", acc, ": ", e$message)
        info <- create_empty_annotation(acc)
      })

      # Always ensure we have something to store
      if (is.null(info)) {
        info <- create_empty_annotation(acc)
      }

      uniprot_data[[acc]] <- info
    }
    else {
      # In offline mode, create empty info
      if (verbose && enable_debug) message("Creating empty annotation for ", acc, " (offline mode)")
      uniprot_data[[acc]] <- create_empty_annotation(acc)
    }

    # Show progress
    if (verbose && (i %% 10 == 0 || i == length(hit_accessions))) {
      message("Processed ", i, " of ", length(hit_accessions), " accessions (",
              cache_hits, " from cache, ", api_calls, " API calls)")
    }
  }

  if (verbose) {
    message("UniProt data retrieval complete: ", cache_hits, " from cache, ",
            api_calls, " from API, ", storage_successes, " stored in database")
  }

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

    if (verbose) message("Processing batch ", ceiling(batch_start/batch_size), " of ",
                         ceiling(nrow(blast_results)/batch_size),
                         " (rows ", batch_start, " to ", batch_end, ")")

    # Start transaction for this batch
    DBI::dbExecute(con, "BEGIN TRANSACTION")

    batch_success <- 0
    batch_ids <- c()

    tryCatch({
      for (i in 1:nrow(current_batch)) {
        blast_result_id <- current_batch$blast_result_id[i]
        hit_accession <- current_batch$hit_accession[i]

        # Enable debug for specific accessions
        enable_debug <- !is.null(debug_accessions) && hit_accession %in% debug_accessions

        # Store annotation using the data we already have
        annotation_id <- store_annotation(
          con,
          blast_result_id,
          uniprot_data[[hit_accession]],
          verbose = enable_debug,
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

  if (verbose) message("Successfully annotated ", successful_annotations, " of ",
                       nrow(blast_results), " BLAST results.")

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

    if (offline_mode) {
      message("\nNOTE: Annotation was performed in offline mode.")
      message("You can re-run the annotation later with UniProt API access to get GO terms and KEGG pathways.")
    }
  }

  # Verify storage if requested
  verification_results <- NULL
  if (verify_storage && length(annotation_ids) > 0) {
    if (verbose) message("Verifying annotation storage...")

    # Check for annotations missing GO terms or KEGG refs
    if (successful_annotations > 0) {
      # Create ID list for query
      id_list <- paste(annotation_ids, collapse = ",")

      # Check for annotations without GO terms
      missing_go_query <- paste0(
        "SELECT a.annotation_id, a.uniprot_accession
         FROM annotations a
         WHERE a.annotation_id IN (", id_list, ")
         AND NOT EXISTS (
           SELECT 1 FROM go_terms g
           WHERE g.annotation_id = a.annotation_id
         )"
      )

      missing_go <- DBI::dbGetQuery(con, missing_go_query)

      # Check for annotations without KEGG references
      missing_kegg_query <- paste0(
        "SELECT a.annotation_id, a.uniprot_accession
         FROM annotations a
         WHERE a.annotation_id IN (", id_list, ")
         AND NOT EXISTS (
           SELECT 1 FROM kegg_references k
           WHERE k.annotation_id = a.annotation_id
         )"
      )

      missing_kegg <- DBI::dbGetQuery(con, missing_kegg_query)

      verification_results <- list(
        verified = length(annotation_ids),
        missing_go = nrow(missing_go),
        missing_kegg = nrow(missing_kegg),
        missing_go_list = if(nrow(missing_go) > 0) missing_go else NULL,
        missing_kegg_list = if(nrow(missing_kegg) > 0) missing_kegg else NULL
      )

      if (verbose) {
        message("Verification complete: ", verification_results$verified, " annotations verified")
        if (verification_results$missing_go > 0) {
          message("- ", verification_results$missing_go, " annotations missing GO terms")
        }
        if (verification_results$missing_kegg > 0) {
          message("- ", verification_results$missing_kegg, " annotations missing KEGG references")
        }
      }
    }
  }

  # Return summary
  result <- list(
    blast_param_id = blast_param_id,
    unique_accessions = length(hit_accessions),
    successful_extractions = length(uniprot_data),
    annotated_results = successful_annotations,
    go_terms = go_count,
    kegg_refs = kegg_count,
    offline_mode = offline_mode,
    cache_used = use_cache,
    cache_stored = store_cache,
    cache_hits = cache_hits,
    api_calls = api_calls
  )

  # Add verification results if available
  if (!is.null(verification_results)) {
    result$verification <- verification_results
  }

  return(result)
}

#' Query the UniProt API for a protein accession
#'
#' @param accession The UniProt accession number
#' @param fields The fields to retrieve from UniProt (for search endpoint)
#' @param debug If TRUE, save the API response to a file for debugging
#'
#' @return A list containing API response with consistent structure
#' @export
query_uniprot_api <- function(accession,
                              fields = "accession,id,gene_names,go,xref_kegg,organism_name,protein_name",
                              debug = FALSE) {
  # Create result structure
  result <- list(
    accession = accession,
    url = NA_character_,
    status_code = NA_integer_,
    content = NA_character_,
    data = NULL,
    error = NA_character_
  )

  # Try direct URL first (more reliable and complete)
  url <- paste0("https://rest.uniprot.org/uniprotkb/", accession, ".json")
  result$url <- url

  if (debug) message("Querying direct URL: ", url)

  # Make HTTP request with proper headers
  headers <- c(
    "User-Agent" = "R/funseqR UniProt Client",
    "Accept" = "application/json"
  )

  # Try direct URL
  response <- httr::GET(
    url,
    httr::add_headers(.headers = headers),
    httr::timeout(30)
  )

  # Get status code
  status <- httr::status_code(response)
  result$status_code <- status

  # If direct URL fails, try search endpoint
  if (status != 200) {
    if (debug) message("Direct URL failed with status ", status, ", trying search endpoint...")

    search_url <- paste0("https://rest.uniprot.org/uniprotkb/search?query=accession:",
                         accession, "&fields=", URLencode(fields), "&format=json")

    response <- httr::GET(
      search_url,
      httr::add_headers(.headers = headers),
      httr::timeout(30)
    )

    status <- httr::status_code(response)
    result$status_code <- status
    result$url <- search_url

    if (debug) message("Search URL status code: ", status)
  }

  # Process response if successful
  if (status == 200) {
    # Get the raw content
    raw_content <- httr::content(response, "raw")

    if (length(raw_content) > 0) {
      if (debug) message("Raw content length: ", length(raw_content), " bytes")

      # Check for gzip compression (magic bytes 0x1f 0x8b)
      is_gzipped <- length(raw_content) > 2 && raw_content[1] == 0x1f && raw_content[2] == 0x8b

      if (is_gzipped) {
        if (debug) message("Content is gzip compressed, decompressing...")

        # Decompress gzipped content
        tryCatch({
          # Save to a temporary file and use gzfile
          temp_gz <- tempfile(fileext = ".gz")
          writeBin(raw_content, temp_gz)

          con <- gzfile(temp_gz, "rb")
          text_lines <- readLines(con, warn = FALSE)
          close(con)
          text_content <- paste(text_lines, collapse = "\n")

          if (debug) {
            message("Decompressed content length: ", nchar(text_content), " chars")
            if (nchar(text_content) > 0) {
              message("First 100 chars: ", substr(text_content, 1, 100))
            }
          }

          # Store decompressed content
          result$content <- text_content

          # Save to debug file if requested
          if (debug) {
            debug_dir <- "uniprot_debug"
            if (!dir.exists(debug_dir)) dir.create(debug_dir)
            json_file <- file.path(debug_dir, paste0("uniprot_", accession, ".json"))
            writeLines(text_content, json_file)
            message("Saved decompressed JSON to ", json_file)
          }

          # Parse the JSON
          if (nchar(text_content) > 0) {
            parsed_data <- jsonlite::fromJSON(text_content, flatten = TRUE)
            result$data <- parsed_data

            if (debug) {
              message("Successfully parsed JSON")
            }
          }
        }, error = function(e) {
          if (debug) message("Error decompressing or parsing content: ", e$message)
          result$error <- paste("Error processing content:", e$message)
        })
      } else {
        # Not compressed, process normally
        tryCatch({
          text_content <- rawToChar(raw_content)
          result$content <- text_content

          if (debug) {
            message("Content length: ", nchar(text_content), " chars")
            if (nchar(text_content) > 0) {
              message("First 100 chars: ", substr(text_content, 1, 100))
            }
          }

          # Parse the JSON if we have content
          if (nchar(text_content) > 0) {
            parsed_data <- jsonlite::fromJSON(text_content, flatten = TRUE)
            result$data <- parsed_data
          }
        }, error = function(e) {
          if (debug) message("Error processing content: ", e$message)
          result$error <- paste("Error processing content:", e$message)
        })
      }
    } else {
      if (debug) message("Empty response content (0 bytes)")
      result$error <- "Empty response content"
    }
  } else {
    if (debug) message("Request failed with status code: ", status)
    result$error <- paste("HTTP request failed with status code:", status)
  }

  return(result)
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


#' Verify annotation storage for processed entries
#'
#' @param con Database connection object
#' @param annotation_ids Vector of annotation IDs to verify
#' @param verbose Print debug information
#'
#' @return A list with verification results
verify_annotations <- function(con, annotation_ids, verbose = FALSE) {
  if (length(annotation_ids) == 0) {
    if (verbose) message("No annotation IDs provided")
    return(list(verified = 0, missing_go = 0, missing_kegg = 0))
  }

  # Create ID list for query
  id_list <- paste(annotation_ids, collapse = ",")

  # Check for annotations without GO terms
  missing_go_query <- paste0(
    "SELECT a.annotation_id, a.uniprot_accession
     FROM annotations a
     WHERE a.annotation_id IN (", id_list, ")
     AND NOT EXISTS (
       SELECT 1 FROM go_terms g
       WHERE g.annotation_id = a.annotation_id
     )"
  )

  missing_go <- DBI::dbGetQuery(con, missing_go_query)

  # Check for annotations without KEGG references
  missing_kegg_query <- paste0(
    "SELECT a.annotation_id, a.uniprot_accession
     FROM annotations a
     WHERE a.annotation_id IN (", id_list, ")
     AND NOT EXISTS (
       SELECT 1 FROM kegg_references k
       WHERE k.annotation_id = a.annotation_id
     )"
  )

  missing_kegg <- DBI::dbGetQuery(con, missing_kegg_query)

  if (verbose) {
    message("Verified ", length(annotation_ids), " annotations")
    if (nrow(missing_go) > 0) {
      message("- ", nrow(missing_go), " annotations missing GO terms")
    }
    if (nrow(missing_kegg) > 0) {
      message("- ", nrow(missing_kegg), " annotations missing KEGG references")
    }
  }

  return(list(
    verified = length(annotation_ids),
    missing_go = nrow(missing_go),
    missing_kegg = nrow(missing_kegg),
    missing_go_list = missing_go,
    missing_kegg_list = missing_kegg
  ))
}
