#' Process BLAST results and retrieve UniProt annotations
#'
#' @param blast_results A list containing BLAST results
#' @param max_hits Maximum number of hits to process per query
#' @param e_value_threshold E-value threshold for filtering BLAST hits
#' @param verbose Logical, whether to print verbose output
#' @param delay Delay between API calls in seconds
#'
#' @importFrom httr GET content status_code
#' @importFrom jsonlite fromJSON
#' @importFrom dplyr %>% group_by filter slice_min ungroup mutate
#' @importFrom progress progress_bar
#'
#' @return A list containing annotated BLAST results, GO terms, KEGG references, and API responses
#' @export
process_blast_results <- function(blast_results, max_hits = 5, e_value_threshold = 1e-10, verbose = FALSE, delay = 1) {
  # Process BLAST results
  blast_df <- blast_results$results %>%
    group_by(qseqid) %>%
    filter(evalue <= e_value_threshold) %>%
    slice_min(order_by = evalue, n = max_hits) %>%
    ungroup()
  
  # Extract unique UniProt accessions
  uniprot_ids <- unique(gsub("sp\\|([^|]+)\\|.*", "\\1", blast_df$sseqid))
  
  # Create a data frame to store API responses and errors
  api_responses <- data.frame(
    accession = character(),
    url = character(),
    status_code = integer(),
    content = character(),
    error = character(),
    stringsAsFactors = FALSE
  )
  
  # Function to make UniProt API call
  get_uniprot_data <- function(accession) {
    base_url <- "https://rest.uniprot.org/uniprotkb/search"
    query <- paste0("accession:", accession)
    fields <- "accession,id,gene_names,go_id,go,xref_kegg"
    
    url <- paste0(base_url, "?query=", URLencode(query), "&fields=", URLencode(fields), "&format=json")
    
    response <- GET(url)
    
    Sys.sleep(delay)  # Add delay between API calls
    
    if (status_code(response) == 200) {
      content <- content(response, "text", encoding = "UTF-8")
      data <- fromJSON(content, flatten = TRUE)
      
      # Log the API response
      api_responses <<- rbind(api_responses, data.frame(
        accession = accession,
        url = url,
        status_code = status_code(response),
        content = content,
        error = NA_character_,
        stringsAsFactors = FALSE
      ))
      
      return(data)
    } else {
      warning(paste("Failed to retrieve data for", accession, "- Status code:", status_code(response)))
      
      # Log the error
      api_responses <<- rbind(api_responses, data.frame(
        accession = accession,
        url = url,
        status_code = status_code(response),
        content = NA_character_,
        error = paste("Failed to retrieve data - Status code:", status_code(response)),
        stringsAsFactors = FALSE
      ))
      
      return(NULL)
    }
  }
  
  # Create progress bar
  pb <- progress_bar$new(
    format = "[:bar] :percent ETA: :eta",
    total = length(uniprot_ids),
    clear = FALSE,
    width = 60
  )
  
  # Function to extract UniProt information
  extract_uniprot_info <- function(uniprot_data) {
    if (!is.data.frame(uniprot_data$results) || nrow(uniprot_data$results) == 0) {
      warning("Unexpected or empty response from UniProt API")
      return(NULL)
    }
    
    result <- uniprot_data$results
    
    # Extract gene names
    gene_names <- c()
    if (!is.null(result$genes[[1]])) {
      gene_names <- c(
        if (!is.null(result$genes[[1]]$geneName.value)) result$genes[[1]]$geneName.value else character(0),
        if (!is.null(result$genes[[1]]$synonyms)) unlist(lapply(result$genes[[1]]$synonyms, function(x) x$value)) else character(0),
        if (!is.null(result$genes[[1]]$orfNames)) unlist(lapply(result$genes[[1]]$orfNames, function(x) x$value)) else character(0)
      )
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
    
    if (!is.null(result$uniProtKBCrossReferences[[1]])) {
      go_refs <- result$uniProtKBCrossReferences[[1]][result$uniProtKBCrossReferences[[1]]$database == "GO", ]
      if (nrow(go_refs) > 0) {
        go_terms <- do.call(rbind, lapply(1:nrow(go_refs), function(i) {
          properties <- setNames(
            go_refs$properties[[i]]$value,
            go_refs$properties[[i]]$key
          )
          data.frame(
            go_id = go_refs$id[i],
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
      stringsAsFactors = FALSE
    )
    
    if (!is.null(result$uniProtKBCrossReferences[[1]])) {
      kegg_refs <- result$uniProtKBCrossReferences[[1]][result$uniProtKBCrossReferences[[1]]$database == "KEGG", ]
      if (nrow(kegg_refs) > 0) {
        kegg_refs <- data.frame(
          kegg_id = kegg_refs$id,
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
  
  # Retrieve UniProt data for all unique accessions
  uniprot_data <- lapply(uniprot_ids, function(id) {
    if (verbose) cat("Retrieving data for", id, "\n")
    pb$tick()
    
    tryCatch({
      data <- get_uniprot_data(id)
      if (is.null(data) || length(data$results) == 0) {
        warning(paste("No data returned for", id))
        return(NULL)
      }
      
      info <- extract_uniprot_info(data)
      if (is.null(info)) {
        warning(paste("Failed to extract info for", id))
        return(NULL)
      }
      
      if (verbose) {
        cat("Extracted info for", id, ":\n")
        cat("  Gene names:", info$gene_names, "\n")
        cat("  GO terms:", nrow(info$go_terms), "\n")
        cat("  KEGG refs:", nrow(info$kegg_refs), "\n")
      }
      
      info
    }, error = function(e) {
      warning(paste("Error retrieving data for", id, ":", e$message))
      
      # Log the error
      api_responses <<- rbind(api_responses, data.frame(
        accession = id,
        url = NA_character_,
        status_code = NA_integer_,
        content = NA_character_,
        error = e$message,
        stringsAsFactors = FALSE
      ))
      
      NULL
    })
  })
  names(uniprot_data) <- uniprot_ids
  
  # Remove NULL entries from uniprot_data
  uniprot_data <- uniprot_data[!sapply(uniprot_data, is.null)]
  
  # Combine BLAST results with UniProt annotations
  annotated_results <- blast_df %>%
    mutate(
      uniprot_accession = gsub("sp\\|([^|]+)\\|.*", "\\1", sseqid),
      gene_names = sapply(uniprot_accession, function(acc) {
        if (!is.null(uniprot_data[[acc]])) uniprot_data[[acc]]$gene_names else NA_character_
      })
    )
  
  # Prepare GO terms
  go_terms <- do.call(rbind, lapply(names(uniprot_data), function(acc) {
    go <- uniprot_data[[acc]]$go_terms
    if (!is.null(go) && nrow(go) > 0) {
      qseqids <- annotated_results$qseqid[annotated_results$uniprot_accession == acc]
      suppressWarnings(
        data.frame(
          qseqid = rep(qseqids, each = nrow(go)),
          uniprot_accession = acc,
          go,
          stringsAsFactors = FALSE,
          row.names = NULL
        )
      )
    } else {
      NULL
    }
  }))
  
  # Prepare KEGG references
  kegg_refs <- do.call(rbind, lapply(names(uniprot_data), function(acc) {
    kegg <- uniprot_data[[acc]]$kegg_refs
    if (!is.null(kegg) && nrow(kegg) > 0) {
      qseqids <- annotated_results$qseqid[annotated_results$uniprot_accession == acc]
      suppressWarnings(
        data.frame(
          qseqid = rep(qseqids, each = nrow(kegg)),
          uniprot_accession = acc,
          kegg,
          stringsAsFactors = FALSE,
          row.names = NULL
        )
      )
    } else {
      NULL
    }
  }))
  
  if (verbose) {
    cat("\nSummary:\n")
    cat("Total UniProt IDs processed:", length(uniprot_ids), "\n")
    cat("Successful annotations:", length(uniprot_data), "\n")
    cat("Failed annotations:", length(uniprot_ids) - length(uniprot_data), "\n")
    if (is.null(go_terms)) {
      cat("No GO terms retrieved\n")
    } else {
      cat("GO terms retrieved for", length(unique(go_terms$uniprot_accession)), "proteins\n")
    }
    if (is.null(kegg_refs)) {
      cat("No KEGG references retrieved\n")
    } else {
      cat("KEGG references retrieved for", length(unique(kegg_refs$uniprot_accession)), "proteins\n")
    }
  }
  
  list(
    annotated_blast = annotated_results,
    go_terms = go_terms,
    kegg_refs = kegg_refs,
    api_responses = api_responses
  )
}
