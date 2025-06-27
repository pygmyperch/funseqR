# EXPORTED

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

#' Process BLAST results and retrieve UniProt annotations
#'
#' This function processes BLAST results from the database and retrieves
#' UniProt annotations for the hit accessions, including GO terms and KEGG pathways.
#' It supports evidence-based filtering of GO terms to control annotation quality.
#'
#' @param con Database connection object.
#' @param blast_param_id The ID of the BLAST parameters.
#' @param max_hits Maximum number of BLAST hits to annotate per genomic locus. Default is 1 (best hit only).
#'   Setting this to 1 ensures each locus gets a single, most confident functional annotation based on 
#'   the lowest e-value, preventing redundant annotations and statistical bias in enrichment analyses.
#'   Use higher values only if multiple hits per locus are specifically needed for your analysis.
#' @param e_value_threshold E-value threshold for filtering BLAST hits. Default is 1e-10.
#' @param batch_size Number of annotations to process in each transaction. Default is 500.
#' @param delay Delay between operations in seconds. Default is 1.
#' @param offline_mode If TRUE, skips UniProt API and uses basic annotations. Default is FALSE.
#' @param use_cache If TRUE, uses cached API responses if available. Default is TRUE.
#' @param store_cache If TRUE, stores API responses in the database. Default is TRUE.
#' @param verify_storage If TRUE, verifies annotation storage after processing. Default is TRUE.
#' @param evidence_keep Character vector of evidence codes to keep for GO annotations. 
#'   Default is c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "TAS", "IC", "IEA", "ISS") which includes
#'   experimental, curated, and computational evidence (balanced approach). Set to NULL to keep all evidence codes.
#'   See Details section for complete evidence code definitions and recommendations.
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE.
#' @param debug_accessions Vector of accessions to debug. Default is NULL.
#'
#' @details 
#' \strong{GO Evidence Code Filtering}
#' 
#' Gene Ontology annotations are supported by evidence codes that indicate the type and 
#' reliability of the evidence. This function supports filtering GO terms based on these 
#' evidence codes to control annotation quality.
#' 
#' \strong{Evidence Code Definitions:}
#' \itemize{
#'   \item \strong{EXP} - Inferred from Experiment: Direct experimental evidence
#'   \item \strong{IDA} - Inferred from Direct Assay: Direct molecular interaction/localization
#'   \item \strong{IPI} - Inferred from Physical Interaction: Protein-protein interactions
#'   \item \strong{IMP} - Inferred from Mutant Phenotype: Loss/gain of function studies
#'   \item \strong{IGI} - Inferred from Genetic Interaction: Epistasis or suppression
#'   \item \strong{IEP} - Inferred from Expression Pattern: Temporal/spatial expression
#'   \item \strong{TAS} - Traceable Author Statement: Curator judgment from literature
#'   \item \strong{IC} - Inferred by Curator: Expert biocuration
#'   \item \strong{IEA} - Inferred from Electronic Annotation: Computational prediction
#'   \item \strong{ISS} - Inferred from Sequence Similarity: Homology-based transfer
#'   \item \strong{ISO} - Inferred from Sequence Orthology: Ortholog-based transfer
#'   \item \strong{ISA} - Inferred from Sequence Alignment: Alignment-based transfer
#'   \item \strong{ISM} - Inferred from Sequence Model: Domain/motif-based prediction
#'   \item \strong{IGC} - Inferred from Genomic Context: Synteny or gene neighborhood
#'   \item \strong{IBA} - Inferred from Biological Aspect of Ancestor: Phylogenetic inference
#'   \item \strong{IBD} - Inferred from Biological Aspect of Descendant: Phylogenetic inference
#'   \item \strong{IKR} - Inferred from Key Residues: Critical amino acid analysis
#'   \item \strong{IRD} - Inferred from Rapid Divergence: Evolutionary rate analysis
#'   \item \strong{RCA} - Inferred from Reviewed Computational Analysis: Manual review of computational prediction
#' }
#' 
#' \strong{Evidence Quality Hierarchy:}
#' \enumerate{
#'   \item \strong{Experimental Evidence} (Highest confidence): EXP, IDA, IPI, IMP, IGI, IEP
#'   \item \strong{Curated Evidence} (High confidence): TAS, IC
#'   \item \strong{Computational Evidence} (Medium confidence): ISS, ISO, ISA, ISM, IGC, IBA, IBD, IKR, IRD, RCA
#'   \item \strong{Electronic Evidence} (Lower confidence): IEA
#' }
#' 
#' \strong{Compound Evidence Codes:}
#' 
#' UniProt often provides compound evidence codes (e.g., "IEA:UniProtKB-SubCell", "IMP:ZFIN") 
#' that include the primary evidence code followed by additional source information. This 
#' function automatically extracts the primary evidence code for filtering.
#' 
#' \strong{Recommendations for Non-Model Organisms:}
#' 
#' For non-model organisms like teleost fish, computational annotations are often the 
#' primary source of functional information:
#' 
#' \itemize{
#'   \item \strong{Balanced (Default - Recommended for most analyses):}
#'     \code{evidence_keep = c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "TAS", "IC", "IEA", "ISS")}
#'     Includes computational predictions and homology transfers. Good balance of coverage and quality.
#'   
#'   \item \strong{Conservative (High confidence only):}
#'     \code{evidence_keep = c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "TAS", "IC")}
#'     Use when you need high-confidence annotations only. May result in fewer GO terms.
#'   
#'   \item \strong{Comprehensive (Maximum coverage):}
#'     \code{evidence_keep = NULL}
#'     Accepts all evidence types. Use when functional coverage is more important than evidence quality.
#'   
#'   \item \strong{Custom filtering:}
#'     You can specify any combination of evidence codes based on your specific requirements.
#' }
#' 
#' \strong{Performance Notes:}
#' 
#' The function uses local caching to minimize API calls. First run will be slower as it 
#' retrieves and caches UniProt data. Subsequent runs will be much faster using cached data.
#'
#' @return A list containing annotation statistics.
#'
#' @examples
#' \dontrun{
#' # Connect to database and get BLAST results
#' con <- connect_funseq_db("analysis.db")
#' blast_results <- perform_blast_db(con, vcf_file_id, db_path, db_name, "blastx")
#' 
#' # Balanced annotation (default - recommended for most analyses)
#' annotation_results <- annotate_blast_results(con, blast_results$blast_param_id)
#' 
#' # Conservative annotation (high-confidence evidence only)
#' annotation_results <- annotate_blast_results(con, blast_results$blast_param_id,
#'   evidence_keep = c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "TAS", "IC"))
#' 
#' # Comprehensive annotation (all evidence types)
#' annotation_results <- annotate_blast_results(con, blast_results$blast_param_id,
#'   evidence_keep = NULL)
#' 
#' # Debug specific problematic accessions
#' annotation_results <- annotate_blast_results(con, blast_results$blast_param_id,
#'   debug_accessions = c("P12345", "Q67890"),
#'   verbose = TRUE)
#' 
#' # Check results
#' str(annotation_results)
#' print(paste("GO terms extracted:", annotation_results$go_terms))
#' }
#'
#' @export
annotate_blast_results <- function(con, blast_param_id, max_hits = 1, e_value_threshold = 1e-10,
                                   batch_size = 500, delay = 1, offline_mode = FALSE,
                                   use_cache = TRUE, store_cache = TRUE, verify_storage = TRUE,
                                   evidence_keep = c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "TAS", "IC", "IEA", "ISS"),
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
    require_uniprot_connection(verbose = verbose)
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
        uniprot_data[[acc]] <- extract_uniprot_info(cache_data, evidence_keep = evidence_keep, debug = enable_debug)
        cache_hits <- cache_hits + 1
        next  # Skip to next accession
      }
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
        info <- extract_uniprot_info(api_result, evidence_keep = evidence_keep, debug = enable_debug)
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

  # Update analysis report if it exists
  tryCatch({
    # Calculate annotation statistics
    total_blast_hits <- length(unique_accessions)
    annotation_rate <- if (total_blast_hits > 0) round((successful_annotations / total_blast_hits) * 100, 1) else 0
    go_per_annotation <- if (successful_annotations > 0) round(go_count / successful_annotations, 1) else 0
    kegg_per_annotation <- if (successful_annotations > 0) round(kegg_count / successful_annotations, 1) else 0

    annotation_message <- paste0(
      "**Functional Annotation Completed**\n\n",
      "- **Unique BLAST hits processed:** ", format(total_blast_hits, big.mark = ","), "\n",
      "- **Successfully annotated:** ", format(successful_annotations, big.mark = ","), " (", annotation_rate, "%)\n",
      "- **Failed annotations:** ", format(failed_annotations, big.mark = ","), "\n",
      "- **GO terms extracted:** ", format(go_count, big.mark = ","), " (avg ", go_per_annotation, " per protein)\n",
      "- **KEGG references:** ", format(kegg_count, big.mark = ","), " (avg ", kegg_per_annotation, " per protein)\n",
      "- **API calls made:** ", format(api_calls, big.mark = ","), "\n",
      "- **Cache hits:** ", format(cache_hits, big.mark = ","), "\n",
      "- **E-value threshold:** ", e_value_threshold, "\n",
      "- **Max hits per query:** ", max_hits
    )

    # TODO: Update this when report functions are fixed
    # update_analysis_report(
    #   con,
    #   section = "annotation",
    #   message = annotation_message,
    #   verbose = FALSE
    # )
  }, error = function(e) {
    # Silently ignore if no report exists
  })

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
#'
#' @return A list containing counts of annotations, GO terms, and KEGG references.
#'
#' @importFrom DBI dbGetQuery
#' @export
count_annotations <- function(con, blast_param_id = NULL) {
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





# INTERNAL

#' Ensure the uniprot_cache table exists
#'
#' @param con A database connection object
#' @param verbose Logical. If TRUE, print progress information
#' @return Invisible NULL
ensure_uniprot_cache_table <- function(con, verbose = FALSE) {
  # Check if table exists
  tables <- DBI::dbListTables(con)
  if (!"uniprot_cache" %in% tables) {
    stop("Database schema is outdated. Please recreate the database with create_funseq_db().")
  }
  if (verbose) message("UniProt cache table verified")
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

    # Ensure the cache table exists
    ensure_uniprot_cache_table(con, verbose = debug)

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

#' connection check for annotate_blast_results function
#'
require_uniprot_connection <- function(verbose = TRUE) {
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
#' @param evidence_keep Character vector of evidence codes to keep for GO annotations.
#'   If NULL, all evidence codes are kept. Default high-confidence codes include
#'   experimental ("EXP", "IDA", "IPI", "IMP", "IGI", "IEP") and curated ("TAS", "IC").
#' @param debug If TRUE, print debugging information
#'
#' @return A list containing extracted information
extract_uniprot_info <- function(uniprot_data, evidence_keep = NULL, debug = FALSE) {
  # Initialize standard result structure
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
    ),
    pfam_domains = data.frame(
      pfam_id = character(0),
      domain_name = character(0),
      match_status = character(0),
      stringsAsFactors = FALSE
    ),
    interpro_families = data.frame(
      interpro_id = character(0),
      family_name = character(0),
      stringsAsFactors = FALSE
    ),
    eggnog_categories = data.frame(
      eggnog_id = character(0),
      taxonomic_scope = character(0),
      stringsAsFactors = FALSE
    )
  )

  # Validate input
  if (is.null(uniprot_data) || !is.list(uniprot_data)) {
    if (debug) message("No valid UniProt data provided")
    return(result)
  }

  # Parse content field (most reliable source)
  entry <- NULL

  if (!is.null(uniprot_data$content) && is.character(uniprot_data$content) && nchar(uniprot_data$content) > 0) {
    if (debug) message("Parsing content field, length: ", nchar(uniprot_data$content))

    tryCatch({
      # CRITICAL FIX: Use simplifyVector = FALSE to preserve nested structures
      parsed <- jsonlite::fromJSON(uniprot_data$content, simplifyVector = FALSE)

      # Find the entry data
      if (!is.null(parsed$results) && is.list(parsed$results) && length(parsed$results) > 0) {
        entry <- parsed$results[[1]]
        if (debug) message("Found entry in results array")
      } else if (!is.null(parsed$primaryAccession)) {
        entry <- parsed
        if (debug) message("Found entry at root level")
      }
    }, error = function(e) {
      if (debug) message("Error parsing JSON content: ", e$message)
    })
  }

  # Try data field if content parsing failed
  if (is.null(entry) && !is.null(uniprot_data$data)) {
    if (debug) message("Trying data field")

    if (is.list(uniprot_data$data)) {
      if (!is.null(uniprot_data$data$results) && is.list(uniprot_data$data$results) &&
          length(uniprot_data$data$results) > 0) {
        entry <- uniprot_data$data$results[[1]]
      } else if (!is.null(uniprot_data$data$primaryAccession)) {
        entry <- uniprot_data$data
      }
    }
  }

  # Return empty result if no entry found
  if (is.null(entry)) {
    if (debug) message("No entry data found in response")
    return(result)
  }

  # Extract basic info
  if (!is.null(entry$primaryAccession)) {
    result$accession <- entry$primaryAccession
    if (debug) message("Found accession: ", result$accession)
  }

  if (!is.null(entry$uniProtkbId)) {
    result$entry_name <- entry$uniProtkbId
    if (debug) message("Found entry name: ", result$entry_name)
  }

  # Extract gene names
  if (!is.null(entry$genes) && is.list(entry$genes) && length(entry$genes) > 0) {
    gene_names <- c()

    for (i in seq_along(entry$genes)) {
      gene <- entry$genes[[i]]
      if (is.list(gene) && !is.null(gene$geneName) &&
          is.list(gene$geneName) && !is.null(gene$geneName$value)) {
        gene_names <- c(gene_names, gene$geneName$value)
        if (debug) message("Found gene name: ", gene$geneName$value)
      }
    }

    if (length(gene_names) > 0) {
      result$gene_names <- paste(gene_names, collapse = ";")
      if (debug) message("Combined gene names: ", result$gene_names)
    }
  }

  # Extract GO terms and KEGG references
  if (!is.null(entry$uniProtKBCrossReferences) && is.list(entry$uniProtKBCrossReferences)) {
    if (debug) message("Processing ", length(entry$uniProtKBCrossReferences), " cross-references")

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

    pfam_domains <- data.frame(
      pfam_id = character(0),
      domain_name = character(0),
      match_status = character(0),
      stringsAsFactors = FALSE
    )

    interpro_families <- data.frame(
      interpro_id = character(0),
      family_name = character(0),
      stringsAsFactors = FALSE
    )

    eggnog_categories <- data.frame(
      eggnog_id = character(0),
      taxonomic_scope = character(0),
      stringsAsFactors = FALSE
    )

    # Process each cross-reference
    for (i in seq_along(entry$uniProtKBCrossReferences)) {
      ref <- entry$uniProtKBCrossReferences[[i]]

      # Process GO terms - FIXED: properly check structure and extract fields
      if (is.list(ref) && !is.null(ref$database) && ref$database == "GO" && !is.null(ref$id)) {
        if (debug) message("Processing GO reference: ", ref$id)

        go_id <- ref$id
        go_term <- NA_character_
        go_evidence <- NA_character_

        # Extract properties - FIXED: properly handle the properties list
        if (!is.null(ref$properties) && is.list(ref$properties)) {
          for (j in seq_along(ref$properties)) {
            prop <- ref$properties[[j]]

            if (is.list(prop) && !is.null(prop$key) && !is.null(prop$value)) {
              if (prop$key == "GoTerm") {
                go_term <- prop$value
                if (debug) message("  Found GO term: ", go_term)
              }
              if (prop$key == "GoEvidenceType") {
                go_evidence <- prop$value
                if (debug) message("  Found GO evidence: ", go_evidence)
              }
            }
          }
        }

        # Add GO term if we have one and it passes evidence filter
        if (!is.na(go_term)) {
          # Apply evidence filter if specified
          if (!is.null(evidence_keep) && !is.na(go_evidence)) {
            # Extract the primary evidence code (before any colon)
            primary_evidence <- strsplit(go_evidence, ":")[[1]][1]
            if (!primary_evidence %in% evidence_keep) {
              if (debug) message("Filtered out GO term with evidence: ", go_evidence, " (primary: ", primary_evidence, ")")
              next
            }
          }
          
          # Extract category (C, F, P) from the term
          go_category <- substr(go_term, 1, 1)

          new_row <- data.frame(
            go_id = go_id,
            go_term = go_term,
            go_category = go_category,
            go_evidence = go_evidence,
            stringsAsFactors = FALSE
          )
          go_terms <- rbind(go_terms, new_row)

          if (debug) message("Added GO term: ", go_id, " - ", go_term, " (evidence: ", go_evidence, ")")
        }
      }

      # Process KEGG references
      if (is.list(ref) && !is.null(ref$database) && ref$database == "KEGG" && !is.null(ref$id)) {
        if (debug) message("Processing KEGG reference: ", ref$id)

        kegg_id <- ref$id
        pathway_name <- NA_character_

        # Extract properties - enhanced to handle multiple property key types
        if (!is.null(ref$properties) && is.list(ref$properties)) {
          if (debug) message("  Processing ", length(ref$properties), " properties")
          for (j in seq_along(ref$properties)) {
            prop <- ref$properties[[j]]

            if (is.list(prop) && !is.null(prop$key) && !is.null(prop$value)) {
              if (debug) message("    Property: ", prop$key, " = ", prop$value)
              # Check for various possible pathway name keys
              if (prop$key %in% c("Description", "PathwayName", "Pathway", "Name", "FullName")) {
                pathway_name <- prop$value
                if (debug) message("  Found pathway name via '", prop$key, "': ", pathway_name)
              }
            }
          }
        } else {
          if (debug) message("  No properties found for KEGG reference")
        }

        # Handle missing pathway names - provide more informative fallback
        if (is.na(pathway_name) || pathway_name == "") {
          # For gene IDs (like dre:563201), extract organism info
          if (grepl("^[a-z]{3}:", kegg_id)) {
            organism_code <- substr(kegg_id, 1, 3)
            organism_map <- c(
              "dre" = "Danio rerio (zebrafish)",
              "tru" = "Takifugu rubripes (fugu)",
              "hsa" = "Homo sapiens (human)",
              "mmu" = "Mus musculus (mouse)",
              "rno" = "Rattus norvegicus (rat)"
            )
            if (organism_code %in% names(organism_map)) {
              pathway_name <- paste0("Gene from ", organism_map[organism_code])
            } else {
              pathway_name <- paste0("Gene from organism: ", organism_code)
            }
            if (debug) message("  Generated descriptive name: ", pathway_name)
          } else {
            pathway_name <- "Unknown pathway"
          }
        }

        # Add KEGG reference
        new_row <- data.frame(
          kegg_id = kegg_id,
          pathway_name = pathway_name,
          stringsAsFactors = FALSE
        )
        kegg_refs <- rbind(kegg_refs, new_row)

        if (debug) message("Added KEGG reference: ", kegg_id, " -> ", pathway_name)
      }

      # Process Pfam domains
      if (is.list(ref) && !is.null(ref$database) && ref$database == "Pfam" && !is.null(ref$id)) {
        if (debug) message("Processing Pfam reference: ", ref$id)

        pfam_id <- ref$id
        domain_name <- NA_character_
        match_status <- NA_character_

        # Extract properties
        if (!is.null(ref$properties) && is.list(ref$properties)) {
          for (j in seq_along(ref$properties)) {
            prop <- ref$properties[[j]]

            if (is.list(prop) && !is.null(prop$key) && !is.null(prop$value)) {
              if (prop$key == "EntryName") {
                domain_name <- prop$value
                if (debug) message("  Found domain name: ", domain_name)
              }
              if (prop$key == "MatchStatus") {
                match_status <- prop$value
                if (debug) message("  Found match status: ", match_status)
              }
            }
          }
        }

        # Add Pfam domain
        new_row <- data.frame(
          pfam_id = pfam_id,
          domain_name = domain_name,
          match_status = match_status,
          stringsAsFactors = FALSE
        )
        pfam_domains <- rbind(pfam_domains, new_row)

        if (debug) message("Added Pfam domain: ", pfam_id, " -> ", domain_name)
      }

      # Process InterPro families
      if (is.list(ref) && !is.null(ref$database) && ref$database == "InterPro" && !is.null(ref$id)) {
        if (debug) message("Processing InterPro reference: ", ref$id)

        interpro_id <- ref$id
        family_name <- NA_character_

        # Extract properties
        if (!is.null(ref$properties) && is.list(ref$properties)) {
          for (j in seq_along(ref$properties)) {
            prop <- ref$properties[[j]]

            if (is.list(prop) && !is.null(prop$key) && !is.null(prop$value)) {
              if (prop$key == "EntryName") {
                family_name <- prop$value
                if (debug) message("  Found family name: ", family_name)
              }
            }
          }
        }

        # Add InterPro family
        new_row <- data.frame(
          interpro_id = interpro_id,
          family_name = family_name,
          stringsAsFactors = FALSE
        )
        interpro_families <- rbind(interpro_families, new_row)

        if (debug) message("Added InterPro family: ", interpro_id, " -> ", family_name)
      }

      # Process eggNOG categories
      if (is.list(ref) && !is.null(ref$database) && ref$database == "eggNOG" && !is.null(ref$id)) {
        if (debug) message("Processing eggNOG reference: ", ref$id)

        eggnog_id <- ref$id
        taxonomic_scope <- NA_character_

        # Extract properties
        if (!is.null(ref$properties) && is.list(ref$properties)) {
          for (j in seq_along(ref$properties)) {
            prop <- ref$properties[[j]]

            if (is.list(prop) && !is.null(prop$key) && !is.null(prop$value)) {
              if (prop$key == "ToxonomicScope") {
                taxonomic_scope <- prop$value
                if (debug) message("  Found taxonomic scope: ", taxonomic_scope)
              }
            }
          }
        }

        # Add eggNOG category
        new_row <- data.frame(
          eggnog_id = eggnog_id,
          taxonomic_scope = taxonomic_scope,
          stringsAsFactors = FALSE
        )
        eggnog_categories <- rbind(eggnog_categories, new_row)

        if (debug) message("Added eggNOG category: ", eggnog_id, " -> ", taxonomic_scope)
      }
    }

    # Update result with extracted data
    result$go_terms <- go_terms
    result$kegg_refs <- kegg_refs
    result$pfam_domains <- pfam_domains
    result$interpro_families <- interpro_families
    result$eggnog_categories <- eggnog_categories

    if (debug) {
      message("Extracted ", nrow(go_terms), " GO terms, ", nrow(kegg_refs), " KEGG references, ",
              nrow(pfam_domains), " Pfam domains, ", nrow(interpro_families), " InterPro families, ",
              "and ", nrow(eggnog_categories), " eggNOG categories")
    }
  } else if (debug) {
    message("No cross-references found or not a list structure")
  }

  return(result)
}

#' Create an empty annotation structure
#'
#' @param accession The accession number
#' @return A list with empty annotation structures
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
    ),
    pfam_domains = data.frame(
      pfam_id = character(0),
      domain_name = character(0),
      match_status = character(0),
      stringsAsFactors = FALSE
    ),
    interpro_families = data.frame(
      interpro_id = character(0),
      family_name = character(0),
      stringsAsFactors = FALSE
    ),
    eggnog_categories = data.frame(
      eggnog_id = character(0),
      taxonomic_scope = character(0),
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

  # Add Pfam domains
  pfam_count <- 0
  if (!is.null(annotation_id) && !is.null(uniprot_info$pfam_domains) &&
      is.data.frame(uniprot_info$pfam_domains) && nrow(uniprot_info$pfam_domains) > 0) {
    if (verbose) message("Adding ", nrow(uniprot_info$pfam_domains), " Pfam domains")

    tryCatch({
      for (i in 1:nrow(uniprot_info$pfam_domains)) {
        # Get values with proper NULL handling
        pfam_id <- uniprot_info$pfam_domains$pfam_id[i]
        if (is.na(pfam_id)) next

        # Check if Pfam domain already exists
        pfam_exists <- DBI::dbGetQuery(
          con,
          "SELECT COUNT(*) AS count FROM pfam_domains WHERE annotation_id = ? AND pfam_id = ?",
          params = list(annotation_id, pfam_id)
        )$count > 0

        if (pfam_exists) {
          if (verbose) message("Pfam domain ", pfam_id, " already exists, skipping")
          next
        }

        # Prepare parameters for insertion
        domain_name <- uniprot_info$pfam_domains$domain_name[i]
        match_status <- uniprot_info$pfam_domains$match_status[i]

        # Handle NULL/NA values
        if (is.na(domain_name)) domain_name <- NULL
        if (is.na(match_status)) match_status <- NULL

        # Build query based on what's available
        fields <- c("annotation_id", "pfam_id")
        values <- c("?", "?")
        insert_params <- list(annotation_id, pfam_id)

        if (!is.null(domain_name)) {
          fields <- c(fields, "domain_name")
          values <- c(values, "?")
          insert_params <- c(insert_params, list(domain_name))
        }

        if (!is.null(match_status)) {
          fields <- c(fields, "match_status")
          values <- c(values, "?")
          insert_params <- c(insert_params, list(match_status))
        }

        # Create and execute the insertion query
        insert_query <- paste0(
          "INSERT INTO pfam_domains (", paste(fields, collapse = ", "), ") ",
          "VALUES (", paste(values, collapse = ", "), ")"
        )

        DBI::dbExecute(con, insert_query, params = insert_params)
        pfam_count <- pfam_count + 1
      }

      if (verbose) message("Added ", pfam_count, " Pfam domains")
    }, error = function(e) {
      warning("Error inserting Pfam domains: ", e$message)
    })
  } else if (verbose) {
    message("No Pfam domains to add")
  }

  # Add InterPro families
  interpro_count <- 0
  if (!is.null(annotation_id) && !is.null(uniprot_info$interpro_families) &&
      is.data.frame(uniprot_info$interpro_families) && nrow(uniprot_info$interpro_families) > 0) {
    if (verbose) message("Adding ", nrow(uniprot_info$interpro_families), " InterPro families")

    tryCatch({
      for (i in 1:nrow(uniprot_info$interpro_families)) {
        # Get values with proper NULL handling
        interpro_id <- uniprot_info$interpro_families$interpro_id[i]
        if (is.na(interpro_id)) next

        # Check if InterPro family already exists
        interpro_exists <- DBI::dbGetQuery(
          con,
          "SELECT COUNT(*) AS count FROM interpro_families WHERE annotation_id = ? AND interpro_id = ?",
          params = list(annotation_id, interpro_id)
        )$count > 0

        if (interpro_exists) {
          if (verbose) message("InterPro family ", interpro_id, " already exists, skipping")
          next
        }

        # Prepare parameters for insertion
        family_name <- uniprot_info$interpro_families$family_name[i]
        if (is.na(family_name)) family_name <- NULL

        # Build query based on what's available
        if (is.null(family_name)) {
          insert_query <- "INSERT INTO interpro_families (annotation_id, interpro_id) VALUES (?, ?)"
          insert_params <- list(annotation_id, interpro_id)
        } else {
          insert_query <- "INSERT INTO interpro_families (annotation_id, interpro_id, family_name) VALUES (?, ?, ?)"
          insert_params <- list(annotation_id, interpro_id, family_name)
        }

        DBI::dbExecute(con, insert_query, params = insert_params)
        interpro_count <- interpro_count + 1
      }

      if (verbose) message("Added ", interpro_count, " InterPro families")
    }, error = function(e) {
      warning("Error inserting InterPro families: ", e$message)
    })
  } else if (verbose) {
    message("No InterPro families to add")
  }

  # Add eggNOG categories
  eggnog_count <- 0
  if (!is.null(annotation_id) && !is.null(uniprot_info$eggnog_categories) &&
      is.data.frame(uniprot_info$eggnog_categories) && nrow(uniprot_info$eggnog_categories) > 0) {
    if (verbose) message("Adding ", nrow(uniprot_info$eggnog_categories), " eggNOG categories")

    tryCatch({
      for (i in 1:nrow(uniprot_info$eggnog_categories)) {
        # Get values with proper NULL handling
        eggnog_id <- uniprot_info$eggnog_categories$eggnog_id[i]
        if (is.na(eggnog_id)) next

        # Check if eggNOG category already exists
        eggnog_exists <- DBI::dbGetQuery(
          con,
          "SELECT COUNT(*) AS count FROM eggnog_categories WHERE annotation_id = ? AND eggnog_id = ?",
          params = list(annotation_id, eggnog_id)
        )$count > 0

        if (eggnog_exists) {
          if (verbose) message("eggNOG category ", eggnog_id, " already exists, skipping")
          next
        }

        # Prepare parameters for insertion
        taxonomic_scope <- uniprot_info$eggnog_categories$taxonomic_scope[i]
        if (is.na(taxonomic_scope)) taxonomic_scope <- NULL

        # Build query based on what's available
        if (is.null(taxonomic_scope)) {
          insert_query <- "INSERT INTO eggnog_categories (annotation_id, eggnog_id) VALUES (?, ?)"
          insert_params <- list(annotation_id, eggnog_id)
        } else {
          insert_query <- "INSERT INTO eggnog_categories (annotation_id, eggnog_id, taxonomic_scope) VALUES (?, ?, ?)"
          insert_params <- list(annotation_id, eggnog_id, taxonomic_scope)
        }

        DBI::dbExecute(con, insert_query, params = insert_params)
        eggnog_count <- eggnog_count + 1
      }

      if (verbose) message("Added ", eggnog_count, " eggNOG categories")
    }, error = function(e) {
      warning("Error inserting eggNOG categories: ", e$message)
    })
  } else if (verbose) {
    message("No eggNOG categories to add")
  }

  return(annotation_id)
}

#' Query the UniProt API for a protein accession
#'
#' @param accession The UniProt accession number
#' @param fields The fields to retrieve from UniProt (for search endpoint)
#' @param debug If TRUE, save the API response to a file for debugging
#'
#' @return A list containing API response with consistent structure
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

#' Process UniProt JSON data into annotations
#'
#' Helper function to extract annotation information from parsed UniProt JSON
#'
#' @param parsed Parsed JSON data from UniProt API
#' @param accession The accession number
#' @param evidence_keep Character vector of evidence codes to keep for GO annotations
#' @param debug Enable debug output
#'
#' @return A structured annotation object
process_uniprot_json <- function(parsed, accession, evidence_keep = NULL, debug = FALSE) {
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
            # Apply evidence filter if specified
            if (!is.null(evidence_keep) && !is.na(go_evidence)) {
              # Extract the primary evidence code (before any colon)
              primary_evidence <- strsplit(go_evidence, ":")[[1]][1]
              if (!primary_evidence %in% evidence_keep) {
                if (debug) message("Filtered out GO term with evidence: ", go_evidence, " (primary: ", primary_evidence, ")")
                next
              }
            }
            
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
                # Check for various possible pathway name keys
                if (prop$key %in% c("Description", "PathwayName", "Pathway", "Name", "FullName")) {
                  pathway_name <- prop$value
                }
              }
            }
          }

          # Handle missing pathway names - provide more informative fallback
          if (is.na(pathway_name) || pathway_name == "") {
            # For gene IDs (like dre:563201), extract organism info
            if (grepl("^[a-z]{3}:", kegg_id)) {
              organism_code <- substr(kegg_id, 1, 3)
              organism_map <- c(
                "dre" = "Danio rerio (zebrafish)",
                "tru" = "Takifugu rubripes (fugu)",
                "hsa" = "Homo sapiens (human)",
                "mmu" = "Mus musculus (mouse)",
                "rno" = "Rattus norvegicus (rat)"
              )
              if (organism_code %in% names(organism_map)) {
                pathway_name <- paste0("Gene from ", organism_map[organism_code])
              } else {
                pathway_name <- paste0("Gene from organism: ", organism_code)
              }
            } else {
              pathway_name <- "Unknown pathway"
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

