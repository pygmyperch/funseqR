#' Annotation storage and retrieval functions for funseqR
#'
#' These functions manage UniProt annotations within the funseqR database.
#'

#' Query the UniProt API for a protein accession
#'
#' This function queries the UniProt API for information about a protein accession.
#'
#' @param accession The UniProt accession number.
#' @param fields The fields to retrieve from UniProt. Default is "accession,id,gene_names,go_id,go,xref_kegg".
#'
#' @return A list containing the API response.
#'
#' @importFrom httr GET content status_code
#' @importFrom jsonlite fromJSON
#' @export
query_uniprot_api <- function(accession, fields = "accession,id,gene_names,go_id,go,xref_kegg") {
  base_url <- "https://rest.uniprot.org/uniprotkb/search"
  query <- paste0("accession:", accession)

  url <- paste0(base_url, "?query=", URLencode(query), "&fields=", URLencode(fields), "&format=json")

  response <- httr::GET(url)

  if (httr::status_code(response) == 200) {
    content <- httr::content(response, "text", encoding = "UTF-8")
    data <- jsonlite::fromJSON(content, flatten = TRUE)

    return(list(
      accession = accession,
      url = url,
      status_code = httr::status_code(response),
      content = content,
      data = data,
      error = NA_character_
    ))
  } else {
    warning(paste("Failed to retrieve data for", accession, "- Status code:", httr::status_code(response)))

    return(list(
      accession = accession,
      url = url,
      status_code = httr::status_code(response),
      content = NA_character_,
      data = NULL,
      error = paste("Failed to retrieve data - Status code:", httr::status_code(response))
    ))
  }
}

#' Test the connection to the UniProt API
#'
#' This function tests if the UniProt API is accessible.
#'
#' @return Logical. TRUE if the connection is successful.
#'
#' @importFrom httr GET status_code
#' @export
test_uniprot_connection <- function() {
  test_url <- "https://rest.uniprot.org"

  tryCatch({
    response <- httr::GET(test_url)
    return(httr::status_code(response) == 200)
  }, error = function(e) {
    return(FALSE)
  })
}

#' Extract UniProt information from API response
#'
#' @param uniprot_data The response from query_uniprot_api
#'
#' @return A list containing extracted information
#'
#' @noRd
extract_uniprot_info <- function(uniprot_data) {
  if (is.null(uniprot_data$data) || is.null(uniprot_data$data$results) || length(uniprot_data$data$results) == 0) {
    warning("No data found for accession ", uniprot_data$accession)
    return(NULL)
  }

  result <- uniprot_data$data$results[[1]]

  # Extract gene names
  gene_names <- c()
  if (!is.null(result$genes)) {
    gene_names <- unlist(lapply(result$genes, function(gene) {
      c(
        if (!is.null(gene$geneName$value)) gene$geneName$value else character(0),
        if (!is.null(gene$synonyms)) sapply(gene$synonyms, function(x) x$value) else character(0),
        if (!is.null(gene$orfNames)) sapply(gene$orfNames, function(x) x$value) else character(0)
      )
    }))
    gene_names <- unique(gene_names[!is.na(gene_names)])
  }
  gene_names <- paste(gene_names, collapse = ";")

  # Extract GO terms
  go_terms <- data.frame(
    go_id = character(0),
    go_term = character(0),
    go_category = character(0),
    go_evidence = character(0),
    stringsAsFactors = FALSE
  )

  if (!is.null(result$uniProtKBCrossReferences)) {
    go_refs <- result$uniProtKBCrossReferences[sapply(result$uniProtKBCrossReferences, function(ref) ref$database == "GO")]
    if (length(go_refs) > 0) {
      go_terms <- do.call(rbind, lapply(go_refs, function(ref) {
        properties <- setNames(
          sapply(ref$properties, function(prop) prop$value),
          sapply(ref$properties, function(prop) prop$key)
        )
        data.frame(
          go_id = ref$id,
          go_term = properties["GoTerm"],
          go_category = substr(properties["GoTerm"], 1, 1),
          go_evidence = properties["GoEvidenceType"],
          stringsAsFactors = FALSE
        )
      }))
    }
  }

  # Extract KEGG references
  kegg_refs <- data.frame(
    kegg_id = character(0),
    pathway_name = NA_character_,
    stringsAsFactors = FALSE
  )

  if (!is.null(result$uniProtKBCrossReferences)) {
    kegg_refs_list <- result$uniProtKBCrossReferences[sapply(result$uniProtKBCrossReferences, function(ref) ref$database == "KEGG")]
    if (length(kegg_refs_list) > 0) {
      kegg_refs <- data.frame(
        kegg_id = sapply(kegg_refs_list, function(ref) ref$id),
        pathway_name = NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }

  list(
    accession = result$primaryAccession,
    entry_name = result$uniProtkbId,
    gene_names = gene_names,
    go_terms = go_terms,
    kegg_refs = kegg_refs
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
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE.
#'
#' @return A list containing annotation statistics.
#'
#' @importFrom DBI dbExecute dbGetQuery
#' @importFrom httr GET status_code
#' @export
annotate_blast_results <- function(con, blast_param_id, max_hits = 5, e_value_threshold = 1e-10,
                                   batch_size = 500, delay = 1, verbose = TRUE) {
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

  # Check UniProt API connection with direct test
  if (verbose) message("Testing connection to UniProt API...")

  # Simple direct test
  uniprot_test <- tryCatch({
    test_url <- "https://rest.uniprot.org"
    response <- httr::GET(test_url)
    status_code <- httr::status_code(response)

    if (status_code != 200) {
      if (verbose) message("Warning: UniProt API returned status code: ", status_code)
    } else {
      if (verbose) message("UniProt API connection successful.")
    }

    # Try a simple query as an additional test
    test_acc <- "P99999"  # Cytochrome C - a well-known protein
    test_query <- paste0("https://rest.uniprot.org/uniprotkb/search?query=accession:", test_acc)
    test_response <- httr::GET(test_query)
    test_status <- httr::status_code(test_response)

    if (test_status != 200) {
      stop("UniProt API query test failed with status code: ", test_status)
    } else {
      if (verbose) message("UniProt API query test successful.")
    }

    TRUE  # Return TRUE if all tests pass
  }, error = function(e) {
    if (verbose) message("Error connecting to UniProt API: ", e$message)
    stop("Cannot connect to UniProt API. Please check your internet connection.")
  })

  # Query UniProt for each accession with simple progress reporting
  if (verbose) message("Querying UniProt API for accessions...")

  uniprot_data <- list()
  for (i in seq_along(hit_accessions)) {
    acc <- hit_accessions[i]

    # Query UniProt API
    result <- query_uniprot_api(acc)
    uniprot_data[[acc]] <- result

    # Show progress
    if (verbose && (i %% 10 == 0 || i == length(hit_accessions))) {
      message("Queried ", i, " of ", length(hit_accessions), " accessions")
    }

    # Add delay between API calls
    Sys.sleep(delay)
  }

  # Extract information from UniProt responses
  if (verbose) message("Extracting information from UniProt responses...")

  uniprot_info <- lapply(uniprot_data, extract_uniprot_info)
  valid_info <- uniprot_info[!sapply(uniprot_info, is.null)]

  if (verbose) message("Successfully extracted information for ", length(valid_info), " of ", length(hit_accessions), " accessions.")

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

        # Skip if no UniProt information available
        if (is.null(uniprot_info[[hit_accession]])) {
          next
        }

        # Store annotation
        annotation_id <- store_annotation(
          con,
          blast_result_id,
          uniprot_info[[hit_accession]],
          verbose = FALSE
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
  }

  # Return summary
  list(
    blast_param_id = blast_param_id,
    unique_accessions = length(hit_accessions),
    successful_extractions = length(valid_info),
    annotated_results = successful_annotations,
    go_terms = go_count,
    kegg_refs = kegg_count
  )
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
