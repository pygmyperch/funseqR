#' Retrieve functional annotations for blast results
#'
#' Given blast results, this function retrieves InterPro domain annotations, GO terms, and KEGG reference IDs from the UniProt REST API.
#' Processes data in batches to manage large queries efficiently.
#' Returns detailed functional annotations and summary information.
#' 
#' @param blast_results list containing blast results with results data frame
#' @param output_file base name for output files (optional)
#' @param batch_size number of accessions to process in each batch
#' @param max_retries number of times to retry failed requests
#' @param retry_delay delay in seconds between retries
#' @param timeout timeout in seconds for each request
#' @return list containing:
#'   \itemize{
#'     \item annotated_blast: new annotations added to original blast results table
#'     \item interpro: detailed InterPro domain annotations
#'     \item go: detailed GO term annotations
#'     \item kegg: detailed KEGG pathway annotations
#'     \item metadata: analysis processing information and other metadata
#'     \item queries: record of each api query and response details
#'   }
#' @importFrom dplyr bind_rows distinct mutate select rename group_by summarize ungroup across left_join as_tibble
#' @importFrom tidyr separate_rows
#' @importFrom stringr str_replace_all str_remove_all str_remove str_extract
#' @importFrom progress progress_bar
#' @importFrom utils write.csv
#' @examples
#' \dontrun{
#' # read blast results from file
#' blast_results <- read_blast_results("blast_output.txt")
#'
#' # retrieve annotations with default parameters
#' annotations <- annotate_blast_results(blast_results)
#'
#' # retrieve annotations with custom parameters and save results
#' annotations <- annotate_blast_results(
#'   blast_results,
#'   output_file = "my_annotations",
#'   batch_size = 1000,
#'   timeout = 600
#' )
#'
#' # access results
#' blast_with_annotations <- annotations$annotated_blast
#' interpro_annotations <- annotations$interpro
#' go_annotations <- annotations$go
#' kegg_annotations <- annotations$kegg
#' }
annotate_blast_results <- function(blast_results, 
                                   output_file = NULL,
                                   batch_size = 500,
                                   max_retries = 3,
                                   retry_delay = 5,
                                   timeout = 300) {
  # get unique uniprot ids from blast results
  uniprot_ids <- unique(gsub("sp\\|([^|]+)\\|.*", "\\1", blast_results$results$sseqid))
  batches <- split(uniprot_ids, ceiling(seq_along(uniprot_ids)/batch_size))
  
  message(sprintf("\nprocessing %d unique blast hits", nrow(blast_results$results)))
  message(sprintf("identified %d unique uniprot accessions", length(uniprot_ids)))
  
  # setup progress bar
  pb <- progress::progress_bar$new(
    format = "Retrieving functional annotations [:bar] :current/:total proteins (:percent) | Elapsed: :elapsed | ETA: :eta",
    total = length(uniprot_ids),
    clear = FALSE,
    width = 100,
    show_after = 0
  )
  
  # initialize results
  all_annotations <- list()
  all_queries <- list()
  proteins_processed <- 0
  
  # redirect output to temporary file
  temp_log <- tempfile()
  sink(temp_log, append = TRUE, split = FALSE)
  
  # process batches
  for (i in seq_along(batches)) {
    batch_success <- FALSE
    attempts <- 0
    current_batch_size <- length(batches[[i]])
    
    while (!batch_success && attempts < max_retries) {
      attempts <- attempts + 1
      
      if (attempts > 1) {
        Sys.sleep(retry_delay)
      }
      
      tryCatch({
        batch_result <- query_uniprot_api(
          batches[[i]], 
          timeout = timeout,
          quiet = TRUE
        )
        
        if (!is.null(batch_result$annotations)) {
          all_annotations[[i]] <- batch_result$annotations
          all_queries[[length(all_queries) + 1]] <- list(
            batch = i,
            attempt = attempts,
            accessions = batches[[i]],
            queries = batch_result$queries
          )
          batch_success <- TRUE
          
          # update progress
          sink()
          proteins_processed <- proteins_processed + current_batch_size
          pb$update(proteins_processed / length(uniprot_ids))
          sink(temp_log, append = TRUE, split = FALSE)
        }
      }, error = function(e) {
        sink()
        message("\nError in batch ", i, ": ", e$message)
        sink(temp_log, append = TRUE, split = FALSE)
        
        all_queries[[length(all_queries) + 1]] <- list(
          batch = i,
          attempt = attempts,
          accessions = batches[[i]],
          error = e$message
        )
      })
    }
    
    if (!batch_success) {
      sink()
      message("\nFailed to process batch ", i, " after ", max_retries, " attempts")
      sink(temp_log, append = TRUE, split = FALSE)
    }
    
    Sys.sleep(1)
  }
  
  sink()
  unlink(temp_log)
  
  # process basic info
  basic_info_dfs <- lapply(all_annotations, function(x) x$basic_info)
  if (length(basic_info_dfs) > 0) {
    basic_info <- bind_rows(basic_info_dfs) %>%
      distinct(uniprot_accession, .keep_all = TRUE) %>%
      mutate(
        gene_names = str_replace_all(gene_names, "\\s+", " "), 
        gene_names = str_replace_all(gene_names, "\\s", "; ")
      )
  }
  
  # process interpro annotations
  interpro_dfs <- lapply(all_annotations, function(x) x$interpro)
  if (length(interpro_dfs) > 0) {
    interpro_annotations <- bind_rows(interpro_dfs) %>%
      select(uniprot_accession, interpro_ids, interpro_descriptions) %>%
      separate_rows(interpro_ids, interpro_descriptions, sep = "; ") %>%
      mutate(
        interpro_descriptions = str_remove_all(interpro_descriptions, '[";.]$')
      ) %>%
      distinct() %>%
      rename(interpro_id = interpro_ids, interpro_description = interpro_descriptions)
    
    interpro_aggregated <- interpro_annotations %>%
      group_by(uniprot_accession) %>%
      summarize(
        interpro_ids = paste(sort(unique(interpro_id)), collapse = "; "),
        .groups = "drop"
      )
  } else {
    interpro_annotations <- NULL
    interpro_aggregated <- NULL
  }
  
  # process go annotations
  go_dfs <- lapply(all_annotations, function(x) x$go)
  if (length(go_dfs) > 0) {
    go_annotations <- bind_rows(go_dfs) %>%
      select(uniprot_accession, go_ids, go_terms) %>%
      separate_rows(go_ids, go_terms, sep = "; ") %>%
      mutate(
        term_id = str_extract(go_terms, "GO:[0-9]+"),
        go_term = str_remove(go_terms, " \\[GO:[0-9]+\\]$")
      ) %>%
      filter(term_id == go_ids) %>%
      select(uniprot_accession, go_id = go_ids, go_term) %>%
      distinct()
    
    go_aggregated <- go_annotations %>%
      group_by(uniprot_accession) %>%
      summarize(
        go_ids = paste(sort(unique(go_id)), collapse = "; "),
        .groups = "drop"
      )
  } else {
    go_annotations <- NULL
    go_aggregated <- NULL
  }
  
  # process kegg annotations
  kegg_dfs <- lapply(all_annotations, function(x) x$kegg)
  if (length(kegg_dfs) > 0) {
    kegg_annotations <- bind_rows(kegg_dfs) %>%
      select(uniprot_accession, kegg_refs) %>%
      separate_rows(kegg_refs, sep = "; ") %>%
      mutate(
        kegg_id = str_remove(kegg_refs, ";$")
      ) %>%
      select(-kegg_refs) %>%
      distinct()
    
    kegg_aggregated <- kegg_annotations %>%
      group_by(uniprot_accession) %>%
      summarize(
        kegg_ids = paste(sort(unique(kegg_id)), collapse = "; "),
        .groups = "drop"
      )
  } else {
    kegg_annotations <- NULL
    kegg_aggregated <- NULL
  }
  
  # combine blast results with annotations
  annotated_blast <- blast_results$results %>%
    mutate(
      uniprot_accession = gsub("sp\\|([^|]+)\\|.*", "\\1", sseqid)
    )
  
  # add basic info and annotations
  if (exists("basic_info")) {
    annotated_blast <- annotated_blast %>%
      left_join(basic_info, by = "uniprot_accession") %>%
      as_tibble()
  }
  
  if (!is.null(interpro_aggregated)) {
    annotated_blast <- annotated_blast %>%
      left_join(interpro_aggregated, by = "uniprot_accession")
  }
  
  if (!is.null(go_aggregated)) {
    annotated_blast <- annotated_blast %>%
      left_join(go_aggregated, by = "uniprot_accession")
  }
  
  if (!is.null(kegg_aggregated)) {
    annotated_blast <- annotated_blast %>%
      left_join(kegg_aggregated, by = "uniprot_accession")
  }
  
  # standardize format
  annotated_blast <- annotated_blast %>%
    as_tibble() %>%
    mutate(across(where(is.character), ~str_replace_all(.x, "\\s*;\\s*", "; ")))
  
  # prepare results
  result <- list(
    annotated_blast = annotated_blast,
    interpro = if(!is.null(interpro_annotations)) as_tibble(interpro_annotations) else NULL,
    go = if(!is.null(go_annotations)) as_tibble(go_annotations) else NULL,
    kegg = if(!is.null(kegg_annotations)) as_tibble(kegg_annotations) else NULL,
    metadata = list(
      date = Sys.time(),
      blast_metadata = blast_results$metadata,
      annotation_parameters = list(
        batch_size = batch_size,
        max_retries = max_retries,
        retry_delay = retry_delay,
        timeout = timeout
      )
    ),
    queries = all_queries
  )
  
  # write output files
  if (!is.null(output_file)) {
    write.csv(result$annotated_blast, 
              file = paste0(output_file, "_annotated_blast.csv"), 
              row.names = FALSE)
    
    if (!is.null(result$interpro)) {
      write.csv(result$interpro, 
                file = paste0(output_file, "_interpro.csv"), 
                row.names = FALSE)
    }
    
    if (!is.null(result$go)) {
      write.csv(result$go, 
                file = paste0(output_file, "_go_terms.csv"), 
                row.names = FALSE)
    }
    
    if (!is.null(result$kegg)) {
      write.csv(result$kegg, 
                file = paste0(output_file, "_kegg_refs.csv"), 
                row.names = FALSE)
    }
    
    # write query log
    query_log <- do.call(rbind, lapply(all_queries, function(q) {
      data.frame(
        batch = q$batch,
        attempt = q$attempt,
        accessions = paste(q$accessions, collapse = ";"),
        url = if (!is.null(q$queries$basic)) q$queries$basic$url else NA,
        query_status = if (!is.null(q$queries$basic)) q$queries$basic$status else NA,
        interpro_status = if (!is.null(q$queries$interpro)) q$queries$interpro$status else NA,
        go_status = if (!is.null(q$queries$go)) q$queries$go$status else NA,
        kegg_status = if (!is.null(q$queries$kegg)) q$queries$kegg$status else NA,
        error = if (!is.null(q$error)) q$error else NA,
        stringsAsFactors = FALSE
      )
    }))
    
    write.csv(query_log,
              file = paste0(output_file, "_query_log.csv"),
              row.names = FALSE)
  }
  
  # print summary
  message("\nResults summary:")
  message(sprintf("Total BLAST hits processed: %d", nrow(blast_results$results)))
  message(sprintf("Unique proteins with annotations: %d", 
                  length(unique(annotated_blast$uniprot_accession))))
  message(sprintf("Unique InterPro annotations: %d", 
                  if (!is.null(result$interpro)) length(unique(result$interpro$interpro_id)) else 0))
  message(sprintf("Unique GO terms: %d", 
                  if (!is.null(result$go)) length(unique(result$go$go_id)) else 0))
  message(sprintf("Unique KEGG references: %d", 
                  if (!is.null(result$kegg)) length(unique(result$kegg$kegg_id)) else 0))
  
  return(result)
}

