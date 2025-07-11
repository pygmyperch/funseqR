#' Compile FunseqR Results into Unified Data Frame
#'
#' Progressive function to compile and enhance analysis results throughout the funseqR workflow.
#' Supports multiple stages: annotations, enrichment, and descriptive analyses, building a 
#' comprehensive results data frame that grows with your analysis.
#'
#' @param con Database connection object
#' @param stage Character. Analysis stage: "annotations", "enrichment", or "descriptive"
#' @param data Data.frame. Existing results to enhance (required for "enrichment" and "descriptive" stages)
#' @param include Character vector. Types of annotations to include: "GO", "KEGG", "Pfam", "InterPro", "eggNOG"
#' @param candidate_loci Character or data.frame. Candidate loci specification (same as process_annotations)
#' @param blast_param_id Integer. Optional. Specific BLAST parameter set to use
#' @param analysis_ids Integer vector. ORA analysis IDs to retrieve (for "enrichment" stage)
#' @param significance_threshold Numeric. FDR threshold for enrichment (default 0.05)
#' @param export_csv Character. Optional. File path to export results as CSV
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return Data frame with results from the specified stage:
#' \itemize{
#'   \item \strong{annotations stage}: Base annotation data (locus_id, annotations, dataset_type, etc.)
#'   \item \strong{enrichment stage}: Adds enrichment columns (enriched, enrichment_fdr, enriched_terms, etc.)
#'   \item \strong{descriptive stage}: Adds descriptive analysis columns (future implementation)
#' }
#'
#' @details
#' This unified function enables transparent, progressive analysis:
#' 
#' \strong{Stage 1 - Annotations:} Extract and process functional annotations
#' \enumerate{
#'   \item Links genomic loci to functional annotations
#'   \item Processes GO, KEGG, and other annotation types
#'   \item Adds candidate/background labeling
#'   \item Creates comprehensive annotation summary
#' }
#' 
#' \strong{Stage 2 - Enrichment:} Add over-representation analysis results
#' \enumerate{
#'   \item Retrieves ORA results from database using analysis_ids
#'   \item Maps enriched terms back to individual loci
#'   \item Adds enrichment status and statistics to each locus
#'   \item Preserves all existing annotation data
#' }
#' 
#' \strong{Enrichment Results Interpretation:}
#' 
#' The enrichment stage marks loci as \code{enriched = TRUE} if they contain genes 
#' appearing in significantly enriched GO terms. This includes both candidate and 
#' background loci, providing two important perspectives:
#' 
#' \itemize{
#'   \item \strong{All enriched loci} (\code{results[results$enriched == TRUE, ]}): 
#'   All genomic positions harboring genes involved in overrepresented biological processes
#'   \item \strong{Candidate enriched loci} (\code{enriched_loci[enriched_loci$dataset_type == "candidate", ]}): 
#'   Specific candidate variants contributing to the enrichment signal
#' }
#' 
#' This distinction is important because candidate loci drove the enrichment signal 
#' (higher proportion vs background), while background loci provide biological context 
#' (these genes/processes exist in the broader dataset too).
#' 
#' \strong{Stage 3 - Descriptive:} Add descriptive analysis results (future)
#' \enumerate{
#'   \item Framework for additional descriptive statistics
#'   \item Functional category summaries
#'   \item Pathway coverage analysis
#' }
#'
#' @examples
#' \dontrun{
#' # Progressive workflow - build results step by step
#' 
#' # Stage 1: Generate base annotation data
#' results <- compile_funseq_results(
#'   con, 
#'   stage = "annotations",
#'   include = c("GO", "KEGG"),
#'   candidate_loci = "candidates.vcf",
#'   export_csv = "step1_annotations.csv"
#' )
#' 
#' # Stage 2: Add enrichment results (after running run_ORA)
#' results <- compile_funseq_results(
#'   con,
#'   stage = "enrichment", 
#'   data = results,
#'   analysis_ids = c(1, 2, 3),
#'   export_csv = "step2_with_enrichment.csv"
#' )
#' 
#' # Extract loci associated with enriched genes
#' enriched_results <- results[results$enriched == TRUE, ]
#' 
#' # Extract only candidate loci associated with enriched genes  
#' candidate_enriched <- enriched_results[enriched_results$dataset_type == "candidate", ]
#' 
#' # View enrichment summary
#' table(results$dataset_type, results$enriched)
#' }
#'
#' @importFrom dplyr group_by summarise first rowwise ungroup mutate left_join
#' @importFrom magrittr %>%
#' @export
compile_funseq_results <- function(con, 
                                  stage = c("annotations", "enrichment", "descriptive"),
                                  data = NULL,
                                  include = c("GO", "KEGG", "Pfam", "InterPro", "eggNOG"),
                                  candidate_loci = NULL,
                                  blast_param_id = NULL,
                                  analysis_ids = NULL,
                                  significance_threshold = 0.05,
                                  export_csv = NULL,
                                  verbose = TRUE) {
  
  # Validate stage parameter
  stage <- match.arg(stage)
  
  if (verbose) message("=== Compiling FunseqR Results: ", stage, " stage ===")
  
  # Route to appropriate stage handler
  if (stage == "annotations") {
    result <- .compile_annotations_stage(con, include, candidate_loci, blast_param_id, verbose)
    
  } else if (stage == "enrichment") {
    if (is.null(data)) {
      stop("'data' parameter is required for enrichment stage. Provide results from annotations stage.")
    }
    if (is.null(analysis_ids)) {
      stop("'analysis_ids' parameter is required for enrichment stage. Provide ORA analysis IDs from run_ORA().")
    }
    result <- .compile_enrichment_stage(con, data, analysis_ids, significance_threshold, verbose)
    
  } else if (stage == "descriptive") {
    if (is.null(data)) {
      stop("'data' parameter is required for descriptive stage. Provide results from previous stages.")
    }
    result <- .compile_descriptive_stage(con, data, verbose)
    
  } else {
    stop("Unknown stage: ", stage)
  }
  
  # Export if requested
  if (!is.null(export_csv)) {
    if (verbose) message("  - Exporting results to: ", export_csv)
    write.csv(result, export_csv, row.names = FALSE)
  }
  
  if (verbose) {
    message("=== ", stage, " stage complete ===")
    message("  - Resulting data frame: ", nrow(result), " rows × ", ncol(result), " columns")
    if (stage == "enrichment" && "enriched" %in% names(result)) {
      enriched_count <- sum(result$enriched, na.rm = TRUE)
      message("  - Loci with enriched terms: ", enriched_count)
    }
  }
  
  return(result)
}

# ==========================
# STAGE HANDLER FUNCTIONS
# ==========================

#' Handle annotations stage
#' @keywords internal
.compile_annotations_stage <- function(con, include, candidate_loci, blast_param_id, verbose) {
  
  if (verbose) message("  - Compiling functional annotations...")
  
  # Validate inputs
  valid_types <- c("GO", "KEGG", "Pfam", "InterPro", "eggNOG")
  include <- match.arg(include, valid_types, several.ok = TRUE)
  
  if (verbose) {
    message("    - Including annotations: ", paste(include, collapse = ", "))
    if (!is.null(blast_param_id)) {
      message("    - Using BLAST parameter set: ", blast_param_id)
    }
  }
  
  # Check database tables
  tables <- DBI::dbListTables(con)
  required_base_tables <- c("vcf_data", "flanking_sequences", "blast_results", "annotations")
  missing_tables <- required_base_tables[!required_base_tables %in% tables]
  
  if (length(missing_tables) > 0) {
    stop("Missing required tables: ", paste(missing_tables, collapse = ", "))
  }
  
  # Build base query for loci with annotations
  base_conditions <- c()
  params <- list()
  
  if (!is.null(blast_param_id)) {
    base_conditions <- c(base_conditions, "bp.blast_param_id = ?")
    params <- append(params, blast_param_id)
  }
  
  base_where <- if (length(base_conditions) > 0) {
    paste0("WHERE ", paste(base_conditions, collapse = " AND "))
  } else {
    ""
  }
  
  if (verbose) message("    - Extracting loci with annotations...")
  
  # Build comprehensive query for all required data
  query <- paste0("
    SELECT DISTINCT
      vd.vcf_id || '_' || vd.chromosome || '_' || vd.position as locus_id,
      vd.chromosome,
      vd.position,
      a.gene_names,
      a.entry_name,
      a.uniprot_accession,
      GROUP_CONCAT(DISTINCT a.annotation_id) as annotation_ids
    FROM vcf_data vd
    JOIN flanking_sequences fs ON vd.vcf_id = fs.vcf_id
    JOIN blast_results br ON fs.flanking_id = br.flanking_id
    JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
    JOIN annotations a ON br.blast_result_id = a.blast_result_id
    ", base_where, "
    GROUP BY vd.vcf_id, vd.chromosome, vd.position, a.gene_names, a.entry_name, a.uniprot_accession
    ORDER BY vd.chromosome, vd.position
  ")
  
  # Execute query with parameters
  result_data <- if (length(params) > 0) {
    DBI::dbGetQuery(con, query, params)
  } else {
    DBI::dbGetQuery(con, query)
  }
  
  if (nrow(result_data) == 0) {
    if (verbose) message("    - No annotated loci found")
    return(data.frame())
  }
  
  if (verbose) message("    - Found ", nrow(result_data), " annotated loci")
  
  # Rename columns to match expected output format
  names(result_data)[names(result_data) == "gene_names"] <- "gene_name"
  names(result_data)[names(result_data) == "entry_name"] <- "protein_name"
  
  # Convert annotation_ids to lists for processing
  result_data$annotation_ids <- lapply(strsplit(result_data$annotation_ids, ","), as.integer)
  
  # Process each annotation type
  if ("GO" %in% include) {
    if (verbose) message("    - Processing GO annotations...")
    result_data <- .process_go_annotations(con, result_data, base_where, params, verbose)
  }
  
  if ("KEGG" %in% include) {
    if (verbose) message("    - Processing KEGG annotations...")
    result_data <- .process_kegg_annotations(con, result_data, base_where, params, verbose)
  }
  
  if ("Pfam" %in% include) {
    if (verbose) message("    - Processing Pfam annotations...")
    result_data <- .process_pfam_annotations(con, result_data, base_where, params, verbose)
  }
  
  if ("InterPro" %in% include) {
    if (verbose) message("    - Processing InterPro annotations...")
    result_data <- .process_interpro_annotations(con, result_data, base_where, params, verbose)
  }
  
  if ("eggNOG" %in% include) {
    if (verbose) message("    - Processing eggNOG annotations...")
    result_data <- .process_eggnog_annotations(con, result_data, base_where, params, verbose)
  }
  
  # Add candidate loci flagging if requested
  if (!is.null(candidate_loci)) {
    if (verbose) message("    - Identifying candidate loci...")
    candidate_locus_ids <- .identify_candidate_loci(con, candidate_loci, verbose)
    result_data$dataset_type <- ifelse(result_data$locus_id %in% candidate_locus_ids, "candidate", "background")
    
    if (verbose) {
      candidate_count <- sum(result_data$dataset_type == "candidate")
      background_count <- sum(result_data$dataset_type == "background")
      message("    - Candidates: ", candidate_count, ", Background: ", background_count)
    }
  } else {
    result_data$dataset_type <- "background"
  }
  
  # Remove the annotation_ids helper column
  result_data$annotation_ids <- NULL
  
  # Summary message
  if (verbose) {
    message("    - Processing completed!")
    message("    - Total loci: ", nrow(result_data))
    
    if ("GO" %in% include) {
      go_count <- sum(!is.na(result_data$go_terms) & result_data$go_terms != "", na.rm = TRUE)
      message("    - Loci with GO annotations: ", go_count)
    }
    if ("KEGG" %in% include) {
      kegg_count <- sum(!is.na(result_data$kegg_pathways) & result_data$kegg_pathways != "", na.rm = TRUE)
      message("    - Loci with KEGG annotations: ", kegg_count)
    }
    if ("Pfam" %in% include) {
      pfam_count <- sum(!is.na(result_data$pfam_domains) & result_data$pfam_domains != "", na.rm = TRUE)
      message("    - Loci with Pfam annotations: ", pfam_count)
    }
    if ("InterPro" %in% include) {
      interpro_count <- sum(!is.na(result_data$interpro_families) & result_data$interpro_families != "", na.rm = TRUE)
      message("    - Loci with InterPro annotations: ", interpro_count)
    }
    if ("eggNOG" %in% include) {
      eggnog_count <- sum(!is.na(result_data$eggnog_categories) & result_data$eggnog_categories != "", na.rm = TRUE)
      message("    - Loci with eggNOG annotations: ", eggnog_count)
    }
  }
  
  return(result_data)
}

#' Handle enrichment stage
#' @keywords internal
.compile_enrichment_stage <- function(con, data, analysis_ids, significance_threshold, verbose) {
  
  if (verbose) message("  - Adding over-representation analysis results...")
  
  # Validate input data
  required_cols <- c("locus_id", "uniprot_accession")
  missing_cols <- required_cols[!required_cols %in% names(data)]
  if (length(missing_cols) > 0) {
    stop("Input data missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Get enrichment data for all analysis IDs
  enrichment_data <- .extract_enrichment_data(con, analysis_ids, significance_threshold, verbose)
  
  if (nrow(enrichment_data) == 0) {
    if (verbose) message("    - No enriched terms found for provided analysis IDs")
    # Add empty enrichment columns
    data$enriched <- FALSE
    data$enrichment_fdr <- NA_real_
    data$enrichment_pvalue <- NA_real_
    data$enriched_terms <- ""
    data$enrichment_analysis_ids <- ""
    return(data)
  }
  
  # Map enriched terms to loci
  enriched_loci_data <- .map_enriched_terms_to_loci(con, data, enrichment_data, verbose)
  
  # Merge enrichment information with original data
  result <- .merge_enrichment_with_data(data, enriched_loci_data, verbose)
  
  return(result)
}

#' Handle descriptive stage (placeholder)
#' @keywords internal
.compile_descriptive_stage <- function(con, data, verbose) {
  
  if (verbose) message("  - Adding descriptive analysis results...")
  
  # Placeholder for future descriptive analysis functionality
  # Could include:
  # - Functional category summaries
  # - Pathway coverage statistics  
  # - Annotation confidence scores
  # - Comparative analysis results
  
  if (verbose) message("    - Descriptive stage not yet implemented")
  
  # For now, just return the data unchanged
  return(data)
}

# ==========================
# ENRICHMENT HELPER FUNCTIONS
# ==========================

#' Extract enrichment data from database
#' @keywords internal
.extract_enrichment_data <- function(con, analysis_ids, significance_threshold, verbose) {
  
  if (verbose) message("    - Retrieving enrichment results from database...")
  
  # First check if analysis_ids exist
  analysis_ids_str <- paste(analysis_ids, collapse = ", ")
  
  check_query <- paste0("
    SELECT analysis_id FROM ora_analyses 
    WHERE analysis_id IN (", analysis_ids_str, ")
  ")
  
  existing_ids <- DBI::dbGetQuery(con, check_query)$analysis_id
  
  if (length(existing_ids) == 0) {
    if (verbose) message("    - No analyses found for IDs: ", analysis_ids_str)
    return(data.frame())
  }
  
  if (length(existing_ids) < length(analysis_ids)) {
    missing_ids <- setdiff(analysis_ids, existing_ids)
    if (verbose) message("    - Warning: Analysis IDs not found: ", paste(missing_ids, collapse = ", "))
  }
  
  # Use only existing IDs
  analysis_ids_str <- paste(existing_ids, collapse = ", ")
  
  # Improved query with data type handling to fix potential data type issues
  enrichment_query <- paste0("
    SELECT DISTINCT
      ora.analysis_id,
      ora.annotation_type,
      ora.term_type,
      res.term_id,
      res.term_name,
      res.p_value,
      res.p_adjusted as fdr,
      res.gene_ids,
      res.fold_enrichment
    FROM ora_analyses ora
    JOIN ora_results res ON ora.analysis_id = res.analysis_id
    WHERE ora.analysis_id IN (", analysis_ids_str, ")
      AND res.p_adjusted IS NOT NULL
      AND res.p_adjusted != ''
      AND CAST(res.p_adjusted AS REAL) <= ?
    ORDER BY CAST(res.p_adjusted AS REAL)
  ")
  
  enrichment_data <- DBI::dbGetQuery(con, enrichment_query, list(significance_threshold))
  
  if (verbose) {
    message("    - Found ", nrow(enrichment_data), " significantly enriched terms")
    if (nrow(enrichment_data) > 0) {
      by_analysis <- table(enrichment_data$analysis_id)
      for (i in names(by_analysis)) {
        message("      - Analysis ID ", i, ": ", by_analysis[i], " terms")
      }
    } else {
      # Debug output when no results found
      debug_query <- paste0("
        SELECT 
          ora.analysis_id,
          COUNT(*) as total_results,
          MIN(CAST(res.p_adjusted AS REAL)) as min_padj,
          MAX(CAST(res.p_adjusted AS REAL)) as max_padj,
          SUM(CASE WHEN CAST(res.p_adjusted AS REAL) <= ? THEN 1 ELSE 0 END) as significant_count
        FROM ora_analyses ora
        JOIN ora_results res ON ora.analysis_id = res.analysis_id
        WHERE ora.analysis_id IN (", analysis_ids_str, ")
          AND res.p_adjusted IS NOT NULL
          AND res.p_adjusted != ''
        GROUP BY ora.analysis_id
      ")
      debug_info <- DBI::dbGetQuery(con, debug_query, list(significance_threshold))
      if (nrow(debug_info) > 0) {
        message("    - Debug: Results by analysis ID:")
        for (i in 1:nrow(debug_info)) {
          row <- debug_info[i, ]
          message("      - Analysis ID ", row$analysis_id, ": ", row$total_results, " total, ", 
                  row$significant_count, " significant (p_adj <= ", significance_threshold, ")")
          message("        Min p_adjusted: ", round(row$min_padj, 6), ", Max p_adjusted: ", round(row$max_padj, 6))
        }
      } else {
        message("    - No results found in ora_results table for these analysis IDs")
      }
    }
  }
  
  return(enrichment_data)
}

#' Map enriched terms to individual loci
#' @keywords internal
.map_enriched_terms_to_loci <- function(con, data, enrichment_data, verbose) {
  
  if (verbose) message("    - Mapping enriched terms to individual loci...")
  
  # Create a data frame to store locus-level enrichment information
  locus_enrichment <- data.frame(
    locus_id = character(0),
    enriched = logical(0),
    enrichment_fdr = numeric(0),
    enrichment_pvalue = numeric(0),
    enriched_terms = character(0),
    enriched_term_ids = character(0),
    enrichment_analysis_ids = character(0),
    stringsAsFactors = FALSE
  )
  
  # Extract all unique genes from enriched terms
  all_enriched_genes <- unique(unlist(strsplit(enrichment_data$gene_ids, "/")))
  all_enriched_genes <- all_enriched_genes[!is.na(all_enriched_genes) & all_enriched_genes != ""]
  
  if (length(all_enriched_genes) == 0) {
    if (verbose) message("    - No genes found in enriched terms")
    return(locus_enrichment)
  }
  
  if (verbose) message("    - Found ", length(all_enriched_genes), " unique genes in enriched terms")
  
  # For each locus in the data, check if it's associated with enriched genes
  for (i in 1:nrow(data)) {
    locus_id <- data$locus_id[i]
    uniprot_acc <- data$uniprot_accession[i]
    
    # Check if this locus has genes that appear in enriched terms
    is_enriched <- FALSE
    matched_terms <- character(0)
    matched_term_ids <- character(0)
    matched_analyses <- character(0)
    best_fdr <- NA_real_
    best_pvalue <- NA_real_
    
    if (!is.na(uniprot_acc) && uniprot_acc != "") {
      # Split uniprot accessions (can be multiple, separated by ;)
      locus_genes <- strsplit(uniprot_acc, ";")[[1]]
      
      # Check if any of the locus genes appear in enriched terms
      for (gene in locus_genes) {
        if (gene %in% all_enriched_genes) {
          is_enriched <- TRUE
          
          # Find which terms contain this gene
          matching_rows <- enrichment_data[grepl(paste0("(^|/)", gene, "(/|$)"), 
                                               enrichment_data$gene_ids), ]
          
          if (nrow(matching_rows) > 0) {
            matched_terms <- c(matched_terms, matching_rows$term_name)
            matched_term_ids <- c(matched_term_ids, matching_rows$term_id)
            matched_analyses <- c(matched_analyses, matching_rows$analysis_id)
            
            # Track best (lowest) FDR and p-value
            if (is.na(best_fdr) || min(matching_rows$fdr) < best_fdr) {
              best_fdr <- min(matching_rows$fdr)
            }
            if (is.na(best_pvalue) || min(matching_rows$p_value) < best_pvalue) {
              best_pvalue <- min(matching_rows$p_value)
            }
          }
        }
      }
    }
    
    # Add to results
    locus_enrichment <- rbind(locus_enrichment, data.frame(
      locus_id = locus_id,
      enriched = is_enriched,
      enrichment_fdr = if (is_enriched) best_fdr else NA_real_,
      enrichment_pvalue = if (is_enriched) best_pvalue else NA_real_,
      enriched_terms = if (is_enriched) paste(unique(matched_terms), collapse = ";") else "",
      enriched_term_ids = if (is_enriched) paste(unique(matched_term_ids), collapse = ";") else "",
      enrichment_analysis_ids = if (is_enriched) paste(unique(matched_analyses), collapse = ";") else "",
      stringsAsFactors = FALSE
    ))
  }
  
  if (verbose) {
    enriched_count <- sum(locus_enrichment$enriched)
    message("    - ", enriched_count, " out of ", nrow(data), " loci are associated with enriched terms")
  }
  
  return(locus_enrichment)
}

#' Merge enrichment data with original data frame
#' @keywords internal
.merge_enrichment_with_data <- function(data, enriched_loci_data, verbose) {
  
  if (verbose) message("    - Merging enrichment results with annotation data...")
  
  # Merge by locus_id
  result <- merge(data, enriched_loci_data, by = "locus_id", all.x = TRUE)
  
  # Fill missing enrichment values for loci not in enriched_loci_data
  result$enriched[is.na(result$enriched)] <- FALSE
  result$enrichment_fdr[is.na(result$enrichment_fdr) & !result$enriched] <- NA_real_
  result$enrichment_pvalue[is.na(result$enrichment_pvalue) & !result$enriched] <- NA_real_
  result$enriched_terms[is.na(result$enriched_terms)] <- ""
  result$enriched_term_ids[is.na(result$enriched_term_ids)] <- ""
  result$enrichment_analysis_ids[is.na(result$enrichment_analysis_ids)] <- ""
  
  # Reorder columns to put enrichment columns at the end
  base_cols <- names(data)
  enrichment_cols <- c("enriched", "enrichment_fdr", "enrichment_pvalue", 
                      "enriched_terms", "enriched_term_ids", "enrichment_analysis_ids")
  result <- result[, c(base_cols, enrichment_cols)]
  
  return(result)
}


#' Process GO annotations for loci
#' @keywords internal
.process_go_annotations <- function(con, result_data, base_where, params, verbose) {
  
  # Get GO terms for each annotation
  go_query <- paste0("
    SELECT 
      a.annotation_id,
      gt.go_id,
      gt.go_term,
      gt.go_category
    FROM annotations a
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
    JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
    JOIN go_terms gt ON a.annotation_id = gt.annotation_id
    ", base_where
  )
  
  go_data <- if (length(params) > 0) {
    DBI::dbGetQuery(con, go_query, params)
  } else {
    DBI::dbGetQuery(con, go_query)
  }
  
  if (nrow(go_data) > 0) {
    # Aggregate GO terms by annotation
    go_summary <- go_data %>%
      group_by(annotation_id) %>%
      summarise(
        go_terms = paste(unique(go_id), collapse = ";"),
        go_names = paste(unique(go_term), collapse = ";"),
        go_categories = paste(unique(go_category), collapse = ";"),
        .groups = "drop"
      )
    
    # Join with result data through annotation IDs
    result_data <- result_data %>%
      rowwise() %>%
      mutate(
        go_terms = {
          matching_go <- go_summary[go_summary$annotation_id %in% annotation_ids, ]
          if (nrow(matching_go) > 0) {
            paste(unique(unlist(strsplit(matching_go$go_terms, ";"))), collapse = ";")
          } else {
            NA_character_
          }
        },
        go_names = {
          matching_go <- go_summary[go_summary$annotation_id %in% annotation_ids, ]
          if (nrow(matching_go) > 0) {
            paste(unique(unlist(strsplit(matching_go$go_names, ";"))), collapse = ";")
          } else {
            NA_character_
          }
        },
        go_categories = {
          matching_go <- go_summary[go_summary$annotation_id %in% annotation_ids, ]
          if (nrow(matching_go) > 0) {
            paste(unique(unlist(strsplit(matching_go$go_categories, ";"))), collapse = ";")
          } else {
            NA_character_
          }
        }
      ) %>%
      ungroup()
    
    if (verbose) {
      go_count <- sum(!is.na(result_data$go_terms) & result_data$go_terms != "", na.rm = TRUE)
      message("    - ", go_count, " loci have GO annotations")
    }
  } else {
    # Add empty GO columns
    result_data$go_terms <- NA_character_
    result_data$go_names <- NA_character_
    result_data$go_categories <- NA_character_
    if (verbose) message("    - No GO annotations found")
  }
  
  return(result_data)
}

#' Process KEGG annotations with enhancements
#' @keywords internal
.process_kegg_annotations <- function(con, result_data, base_where, params, verbose) {
  
  # Get KEGG references for each annotation
  kegg_query <- paste0("
    SELECT 
      a.annotation_id,
      kr.kegg_id,
      kr.pathway_name
    FROM annotations a
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
    JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
    JOIN kegg_references kr ON a.annotation_id = kr.annotation_id
    ", base_where
  )
  
  kegg_data <- if (length(params) > 0) {
    DBI::dbGetQuery(con, kegg_query, params)
  } else {
    DBI::dbGetQuery(con, kegg_query)
  }
  
  if (nrow(kegg_data) > 0) {
    if (verbose) message("    - Enhancing KEGG pathway names with KEGGREST...")
    
    # Update pathway names using KEGGREST
    kegg_data <- .enhance_kegg_pathway_names(kegg_data, verbose)
    
    if (verbose) message("    - Adding KEGG BRITE classifications...")
    
    # Add KEGG BRITE classifications
    kegg_data <- .add_kegg_brite_classifications(kegg_data, verbose)
    
    if (verbose) message("    - Assigning functional modules...")
    
    # Add functional module assignments
    kegg_data <- .add_kegg_functional_modules(kegg_data, verbose)
    
    # Aggregate KEGG information by annotation
    kegg_summary <- kegg_data %>%
      group_by(annotation_id) %>%
      summarise(
        kegg_pathways = paste(unique(kegg_id), collapse = ";"),
        kegg_pathway_names = paste(unique(pathway_name), collapse = ";"),
        kegg_brite_categories = paste(unique(na.omit(brite_category)), collapse = ";"),
        kegg_modules = paste(unique(na.omit(functional_module)), collapse = ";"),
        .groups = "drop"
      )
    
    # Join with result data through annotation IDs
    result_data <- result_data %>%
      rowwise() %>%
      mutate(
        kegg_pathways = {
          matching_kegg <- kegg_summary[kegg_summary$annotation_id %in% annotation_ids, ]
          if (nrow(matching_kegg) > 0) {
            paste(unique(unlist(strsplit(matching_kegg$kegg_pathways, ";"))), collapse = ";")
          } else {
            NA_character_
          }
        },
        kegg_pathway_names = {
          matching_kegg <- kegg_summary[kegg_summary$annotation_id %in% annotation_ids, ]
          if (nrow(matching_kegg) > 0) {
            paste(unique(unlist(strsplit(matching_kegg$kegg_pathway_names, ";"))), collapse = ";")
          } else {
            NA_character_
          }
        },
        kegg_brite_categories = {
          matching_kegg <- kegg_summary[kegg_summary$annotation_id %in% annotation_ids, ]
          if (nrow(matching_kegg) > 0) {
            cats <- unique(unlist(strsplit(matching_kegg$kegg_brite_categories, ";")))
            cats <- cats[!is.na(cats) & cats != ""]
            if (length(cats) > 0) paste(cats, collapse = ";") else NA_character_
          } else {
            NA_character_
          }
        },
        kegg_modules = {
          matching_kegg <- kegg_summary[kegg_summary$annotation_id %in% annotation_ids, ]
          if (nrow(matching_kegg) > 0) {
            mods <- unique(unlist(strsplit(matching_kegg$kegg_modules, ";")))
            mods <- mods[!is.na(mods) & mods != ""]
            if (length(mods) > 0) paste(mods, collapse = ";") else NA_character_
          } else {
            NA_character_
          }
        }
      ) %>%
      ungroup()
    
    if (verbose) {
      kegg_count <- sum(!is.na(result_data$kegg_pathways) & result_data$kegg_pathways != "", na.rm = TRUE)
      message("    - ", kegg_count, " loci have KEGG annotations")
    }
  } else {
    # Add empty KEGG columns
    result_data$kegg_pathways <- NA_character_
    result_data$kegg_pathway_names <- NA_character_
    result_data$kegg_brite_categories <- NA_character_
    result_data$kegg_modules <- NA_character_
    if (verbose) message("    - No KEGG annotations found")
  }
  
  return(result_data)
}

#' Enhance KEGG pathway names using KEGGREST
#' @keywords internal
.enhance_kegg_pathway_names <- function(kegg_data, verbose) {
  
  # Check if KEGGREST is available
  if (!requireNamespace("KEGGREST", quietly = TRUE)) {
    if (verbose) message("      - KEGGREST not available, using existing pathway names")
    return(kegg_data)
  }
  
  # Find pathways needing name enhancement
  needs_update <- is.na(kegg_data$pathway_name) | 
                  kegg_data$pathway_name == "" | 
                  kegg_data$pathway_name == "-" |
                  grepl("^Gene:", kegg_data$pathway_name)
  
  if (sum(needs_update) == 0) {
    if (verbose) message("      - All pathway names are already descriptive")
    return(kegg_data)
  }
  
  if (verbose) message("      - Updating ", sum(needs_update), " pathway names...")
  
  # Process gene IDs to get pathway information
  gene_ids <- unique(kegg_data$kegg_id[needs_update & grepl("^[a-z]{3}:\\d+", kegg_data$kegg_id)])
  
  for (gene_id in gene_ids[1:min(20, length(gene_ids))]) {  # Limit API calls
    tryCatch({
      gene_info <- KEGGREST::keggGet(gene_id)
      if (length(gene_info) > 0 && !is.null(gene_info[[1]]$PATHWAY)) {
        pathways <- gene_info[[1]]$PATHWAY
        pathway_descriptions <- paste(names(pathways), pathways, sep = ": ")
        
        # Update all matching records
        kegg_data[kegg_data$kegg_id == gene_id & needs_update, "pathway_name"] <- 
          paste(pathway_descriptions, collapse = "; ")
      }
      Sys.sleep(0.2)  # Rate limiting
    }, error = function(e) {
      if (verbose) message("        - Failed to get pathways for ", gene_id)
    })
  }
  
  return(kegg_data)
}

#' Add KEGG BRITE classifications using KEGGREST
#' @keywords internal
.add_kegg_brite_classifications <- function(kegg_data, verbose) {
  
  # For now, use a simplified mapping based on pathway IDs
  # This could be enhanced with actual KEGG BRITE API calls
  
  kegg_data$brite_category <- NA_character_
  
  # Extract pathway IDs from pathway names or kegg_ids
  pathway_patterns <- list(
    "Metabolism" = c("00010", "00020", "00030", "00040", "00051", "00052", "00053", 
                     "00500", "00520", "00620", "00630", "00640", "00650", "00660"),
    "Genetic Information Processing" = c("03010", "03013", "03015", "03018", "03020", 
                                        "03022", "03030", "03040", "03050", "03060"),
    "Environmental Information Processing" = c("02010", "02020", "02024", "02025", 
                                              "02026", "02030", "02040", "02060"),
    "Cellular Processes" = c("04110", "04111", "04112", "04113", "04114", "04120", 
                            "04130", "04136", "04137", "04140", "04141", "04142"),
    "Organismal Systems" = c("04010", "04014", "04015", "04020", "04022", "04024", 
                            "04080", "04142", "04144", "04145", "04146", "04150")
  )
  
  for (category in names(pathway_patterns)) {
    patterns <- pathway_patterns[[category]]
    for (pattern in patterns) {
      matches <- grepl(pattern, kegg_data$kegg_id) | grepl(pattern, kegg_data$pathway_name)
      kegg_data$brite_category[matches & is.na(kegg_data$brite_category)] <- category
    }
  }
  
  return(kegg_data)
}

#' Add functional module assignments
#' @keywords internal  
.add_kegg_functional_modules <- function(kegg_data, verbose) {
  
  # Use simplified functional module classification
  module_patterns <- list(
    "Carbohydrate Metabolism" = c("glycolysis", "gluconeogenesis", "pentose phosphate", 
                                 "citrate cycle", "pyruvate", "starch", "sucrose"),
    "Energy Metabolism" = c("oxidative phosphorylation", "photosynthesis", 
                           "carbon fixation", "nitrogen metabolism"),
    "Lipid Metabolism" = c("fatty acid", "steroid", "glycerolipid", "bile acid"),
    "Amino Acid Metabolism" = c("alanine", "arginine", "aspartate", "glycine", 
                               "leucine", "lysine", "phenylalanine"),
    "Signal Transduction" = c("MAPK", "calcium", "mTOR", "Wnt", "Notch", "insulin"),
    "Cell Cycle" = c("cell cycle", "DNA replication", "p53"),
    "Immune System" = c("complement", "toll-like", "NOD-like", "chemokine")
  )
  
  kegg_data$functional_module <- NA_character_
  
  for (module in names(module_patterns)) {
    patterns <- module_patterns[[module]]
    for (pattern in patterns) {
      matches <- grepl(pattern, kegg_data$pathway_name, ignore.case = TRUE)
      kegg_data$functional_module[matches & is.na(kegg_data$functional_module)] <- module
    }
  }
  
  return(kegg_data)
}

#' Process Pfam annotations for loci
#' @keywords internal
.process_pfam_annotations <- function(con, result_data, base_where, params, verbose) {
  
  # Get Pfam domains for each annotation
  pfam_query <- paste0("
    SELECT 
      a.annotation_id,
      pd.pfam_id,
      pd.domain_name,
      pd.match_status
    FROM annotations a
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
    JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
    JOIN pfam_domains pd ON a.annotation_id = pd.annotation_id
    ", base_where
  )
  
  pfam_data <- if (length(params) > 0) {
    DBI::dbGetQuery(con, pfam_query, params)
  } else {
    DBI::dbGetQuery(con, pfam_query)
  }
  
  if (nrow(pfam_data) > 0) {
    # Aggregate Pfam domains by annotation
    pfam_summary <- pfam_data %>%
      group_by(annotation_id) %>%
      summarise(
        pfam_domains = paste(unique(pfam_id), collapse = ";"),
        pfam_domain_names = paste(unique(na.omit(domain_name)), collapse = ";"),
        .groups = "drop"
      )
    
    # Join with result data through annotation IDs
    result_data <- result_data %>%
      rowwise() %>%
      mutate(
        pfam_domains = {
          matching_pfam <- pfam_summary[pfam_summary$annotation_id %in% annotation_ids, ]
          if (nrow(matching_pfam) > 0) {
            paste(unique(unlist(strsplit(matching_pfam$pfam_domains, ";"))), collapse = ";")
          } else {
            NA_character_
          }
        },
        pfam_domain_names = {
          matching_pfam <- pfam_summary[pfam_summary$annotation_id %in% annotation_ids, ]
          if (nrow(matching_pfam) > 0) {
            names <- unique(unlist(strsplit(matching_pfam$pfam_domain_names, ";")))
            names <- names[!is.na(names) & names != ""]
            if (length(names) > 0) paste(names, collapse = ";") else NA_character_
          } else {
            NA_character_
          }
        }
      ) %>%
      ungroup()
    
    if (verbose) {
      pfam_count <- sum(!is.na(result_data$pfam_domains) & result_data$pfam_domains != "", na.rm = TRUE)
      message("    - ", pfam_count, " loci have Pfam annotations")
    }
  } else {
    # Add empty Pfam columns
    result_data$pfam_domains <- NA_character_
    result_data$pfam_domain_names <- NA_character_
    if (verbose) message("    - No Pfam annotations found")
  }
  
  return(result_data)
}

#' Process InterPro annotations for loci
#' @keywords internal
.process_interpro_annotations <- function(con, result_data, base_where, params, verbose) {
  
  # Get InterPro families for each annotation
  interpro_query <- paste0("
    SELECT 
      a.annotation_id,
      if.interpro_id,
      if.family_name
    FROM annotations a
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
    JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
    JOIN interpro_families if ON a.annotation_id = if.annotation_id
    ", base_where
  )
  
  interpro_data <- if (length(params) > 0) {
    DBI::dbGetQuery(con, interpro_query, params)
  } else {
    DBI::dbGetQuery(con, interpro_query)
  }
  
  if (nrow(interpro_data) > 0) {
    # Aggregate InterPro families by annotation
    interpro_summary <- interpro_data %>%
      group_by(annotation_id) %>%
      summarise(
        interpro_families = paste(unique(interpro_id), collapse = ";"),
        interpro_family_names = paste(unique(na.omit(family_name)), collapse = ";"),
        .groups = "drop"
      )
    
    # Join with result data through annotation IDs
    result_data <- result_data %>%
      rowwise() %>%
      mutate(
        interpro_families = {
          matching_interpro <- interpro_summary[interpro_summary$annotation_id %in% annotation_ids, ]
          if (nrow(matching_interpro) > 0) {
            paste(unique(unlist(strsplit(matching_interpro$interpro_families, ";"))), collapse = ";")
          } else {
            NA_character_
          }
        },
        interpro_family_names = {
          matching_interpro <- interpro_summary[interpro_summary$annotation_id %in% annotation_ids, ]
          if (nrow(matching_interpro) > 0) {
            names <- unique(unlist(strsplit(matching_interpro$interpro_family_names, ";")))
            names <- names[!is.na(names) & names != ""]
            if (length(names) > 0) paste(names, collapse = ";") else NA_character_
          } else {
            NA_character_
          }
        }
      ) %>%
      ungroup()
    
    if (verbose) {
      interpro_count <- sum(!is.na(result_data$interpro_families) & result_data$interpro_families != "", na.rm = TRUE)
      message("    - ", interpro_count, " loci have InterPro annotations")
    }
  } else {
    # Add empty InterPro columns
    result_data$interpro_families <- NA_character_
    result_data$interpro_family_names <- NA_character_
    if (verbose) message("    - No InterPro annotations found")
  }
  
  return(result_data)
}

#' Identify candidate loci from various input formats
#' @keywords internal
.identify_candidate_loci <- function(con, candidate_input, verbose = FALSE) {
  
  if (is.null(candidate_input)) {
    return(character(0))
  }
  
  candidate_loci <- character(0)
  
  # Handle different input types
  if (is.character(candidate_input) && length(candidate_input) == 1) {
    # Assume it's a file path
    if (file.exists(candidate_input)) {
      # Check file extension
      if (grepl("\\.vcf$", candidate_input, ignore.case = TRUE)) {
        # VCF file - get file_id and extract loci
        if (verbose) message("    - Processing VCF file: ", basename(candidate_input))
        
        # Check if file is already registered
        file_check <- DBI::dbGetQuery(con, 
          "SELECT file_id FROM input_files WHERE file_name = ? AND file_type = 'vcf'",
          list(basename(candidate_input)))
        
        if (nrow(file_check) > 0) {
          file_id <- file_check$file_id[1]
          
          # Extract coordinates for this file
          coords_query <- "SELECT DISTINCT chromosome, position FROM vcf_data WHERE file_id = ?"
          coords_result <- DBI::dbGetQuery(con, coords_query, list(file_id))
          
          if (verbose) message("    - Found ", nrow(coords_result), " candidate loci from VCF")
          
          # Use coordinate-based matching to find corresponding loci in database
          candidate_loci <- .match_coordinates_to_loci(con, coords_result, verbose)
          
        } else {
          warning("VCF file '", basename(candidate_input), "' not found in database. Use import_vcf() first.")
        }
        
      } else if (grepl("\\.(bed|txt|csv)$", candidate_input, ignore.case = TRUE)) {
        # BED/coordinate file
        if (verbose) message("    - Processing coordinate file: ", basename(candidate_input))
        
        coords <- read.table(candidate_input, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
        
        if (ncol(coords) >= 3) {
          # BED format: chr, start, end
          coords_df <- data.frame(
            chromosome = coords[, 1],
            position = coords[, 2],  # Use start position
            stringsAsFactors = FALSE
          )
        } else if (ncol(coords) >= 2) {
          # Simple format: chr, position
          coords_df <- data.frame(
            chromosome = coords[, 1],
            position = coords[, 2],
            stringsAsFactors = FALSE
          )
        } else {
          stop("Coordinate file must have at least 2 columns (chromosome, position)")
        }
        
        candidate_loci <- .match_coordinates_to_loci(con, coords_df, verbose)
      }
    } else {
      stop("File not found: ", candidate_input)
    }
    
  } else if (is.data.frame(candidate_input)) {
    # Data frame with coordinates
    if (verbose) message("    - Processing coordinate data frame")
    
    required_cols <- c("chromosome", "position")
    if (!all(required_cols %in% names(candidate_input))) {
      stop("Data frame must contain 'chromosome' and 'position' columns")
    }
    
    candidate_loci <- .match_coordinates_to_loci(con, candidate_input, verbose)
    
  } else if (is.character(candidate_input) && length(candidate_input) > 1) {
    # Vector of locus_ids already in correct format
    if (verbose) message("    - Using provided locus_id vector")
    candidate_loci <- candidate_input
    
  } else {
    stop("Unsupported candidate_input format. Use VCF file path, coordinate data frame, or locus_id vector.")
  }
  
  return(unique(candidate_loci))
}

#' Match coordinate data frame to database loci
#' @keywords internal
.match_coordinates_to_loci <- function(con, coords_df, verbose = FALSE) {
  
  if (nrow(coords_df) == 0) {
    return(character(0))
  }
  
  # Create a temporary table for efficient matching
  temp_table <- paste0("temp_candidates_", sample(1000:9999, 1))
  
  # Create temporary table
  create_temp_sql <- paste0("CREATE TEMP TABLE ", temp_table, " (
    chromosome TEXT,
    position INTEGER
  )")
  
  DBI::dbExecute(con, create_temp_sql)
  
  # Insert coordinates
  insert_sql <- paste0("INSERT INTO ", temp_table, " (chromosome, position) VALUES (?, ?)")
  
  for (i in 1:nrow(coords_df)) {
    DBI::dbExecute(con, insert_sql, list(coords_df$chromosome[i], coords_df$position[i]))
  }
  
  # Match to existing loci
  match_query <- paste0("
    SELECT DISTINCT vd.vcf_id || '_' || vd.chromosome || '_' || vd.position as locus_id
    FROM vcf_data vd
    JOIN ", temp_table, " tc ON vd.chromosome = tc.chromosome AND vd.position = tc.position
  ")
  
  loci_result <- DBI::dbGetQuery(con, match_query)
  
  # Clean up temporary table
  DBI::dbExecute(con, paste0("DROP TABLE ", temp_table))
  
  if (verbose) message("    - Matched ", nrow(loci_result), " coordinates to database loci")
  
  return(loci_result$locus_id)
}

#' Process eggNOG annotations for loci
#' @keywords internal
.process_eggnog_annotations <- function(con, result_data, base_where, params, verbose) {
  
  # Get eggNOG categories for each annotation
  eggnog_query <- paste0("
    SELECT 
      a.annotation_id,
      ec.eggnog_id,
      ec.taxonomic_scope
    FROM annotations a
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
    JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
    JOIN eggnog_categories ec ON a.annotation_id = ec.annotation_id
    ", base_where
  )
  
  eggnog_data <- if (length(params) > 0) {
    DBI::dbGetQuery(con, eggnog_query, params)
  } else {
    DBI::dbGetQuery(con, eggnog_query)
  }
  
  if (nrow(eggnog_data) > 0) {
    # Aggregate eggNOG categories by annotation
    eggnog_summary <- eggnog_data %>%
      group_by(annotation_id) %>%
      summarise(
        eggnog_categories = paste(unique(eggnog_id), collapse = ";"),
        eggnog_taxonomic_scopes = paste(unique(na.omit(taxonomic_scope)), collapse = ";"),
        .groups = "drop"
      )
    
    # Join with result data through annotation IDs
    result_data <- result_data %>%
      rowwise() %>%
      mutate(
        eggnog_categories = {
          matching_eggnog <- eggnog_summary[eggnog_summary$annotation_id %in% annotation_ids, ]
          if (nrow(matching_eggnog) > 0) {
            paste(unique(unlist(strsplit(matching_eggnog$eggnog_categories, ";"))), collapse = ";")
          } else {
            NA_character_
          }
        },
        eggnog_taxonomic_scopes = {
          matching_eggnog <- eggnog_summary[eggnog_summary$annotation_id %in% annotation_ids, ]
          if (nrow(matching_eggnog) > 0) {
            scopes <- unique(unlist(strsplit(matching_eggnog$eggnog_taxonomic_scopes, ";")))
            scopes <- scopes[!is.na(scopes) & scopes != ""]
            if (length(scopes) > 0) paste(scopes, collapse = ";") else NA_character_
          } else {
            NA_character_
          }
        }
      ) %>%
      ungroup()
    
    if (verbose) {
      eggnog_count <- sum(!is.na(result_data$eggnog_categories) & result_data$eggnog_categories != "", na.rm = TRUE)
      message("    - ", eggnog_count, " loci have eggNOG annotations")
    }
  } else {
    # Add empty eggNOG columns
    result_data$eggnog_categories <- NA_character_
    result_data$eggnog_taxonomic_scopes <- NA_character_
    if (verbose) message("    - No eggNOG annotations found")
  }
  
  return(result_data)
}

#' Get comprehensive annotations for specific loci
#'
#' Retrieve all available annotation information for specified genomic loci, including
#' functional annotations, enrichment status, and pathway information. Useful for detailed
#' investigation of specific positions identified in analyses.
#'
#' @param con Database connection object
#' @param loci Character vector of loci in "chromosome:position" format, or data.frame with chromosome and position columns. If NULL, bed_file must be provided.
#' @param bed_file Character. Path to BED file containing loci to query. If provided, overrides loci parameter.
#' @param include_enrichment Logical. Whether to include enrichment analysis results. Default is TRUE.
#' @param verbose Logical. Print progress information. Default is TRUE.
#'
#' @return Data.frame with comprehensive annotation information for each locus.
#'         Always includes all possible columns from compile_funseq_results() output,
#'         with NA values for missing data.
#'
#' @details
#' This function provides a convenient way to investigate specific loci in detail.
#' It returns a standardized data.frame with the same column structure as
#' compile_funseq_results(), ensuring consistency across the funseqR workflow.
#'
#' \\strong{Input Formats:}
#' \\itemize{
#'   \\item Single locus: "LG4:3814415"
#'   \\item Multiple loci: c("LG4:3814415", "LG8:22215908", "LG10:28936085")
#'   \\item Data.frame: data.frame(chromosome = c("LG4", "LG8"), position = c(3814415, 22215908))
#'   \\item BED file: Any standard BED format file with chromosome and position information
#' }
#'
#' \\strong{Output Structure:}
#' The function always returns the complete column set, populating NA for missing data:
#' locus_id, chromosome, position, gene_name, protein_name, uniprot_accession,
#' go_terms, go_names, go_categories, kegg_pathways, kegg_pathway_names,
#' kegg_brite_categories, kegg_modules, pfam_domains, pfam_domain_names,
#' interpro_families, interpro_family_names, eggnog_categories,
#' eggnog_taxonomic_scopes, dataset_type, enriched, enrichment_fdr,
#' enrichment_pvalue, enriched_terms, enriched_term_ids, enrichment_analysis_ids.
#'
#' @examples
#' \\dontrun{
#' con <- connect_funseq_db("analysis.db")
#'
#' # Single locus investigation
#' peak_info <- get_locus_annotations(con, "LG4:3814415")
#'
#' # Multiple loci analysis
#' top_candidates <- c("LG4:3814415", "LG8:22215908", "LG10:28936085")
#' candidate_info <- get_locus_annotations(con, top_candidates)
#'
#' # BED file analysis
#' region_info <- get_locus_annotations(con, bed_file = "interesting_regions.bed")
#'
#' # Data.frame input
#' loci_df <- data.frame(chromosome = c("LG4", "LG8"), position = c(3814415, 22215908))
#' batch_info <- get_locus_annotations(con, loci_df)
#'
#' # View specific columns
#' print(candidate_info[, c("chromosome", "position", "gene_name", "enriched", "enriched_terms")])
#' }
#'
#' @export
get_locus_annotations <- function(con, loci = NULL, bed_file = NULL, include_enrichment = TRUE, verbose = TRUE) {
  
  if (is.null(loci) && is.null(bed_file)) {
    stop("Either 'loci' or 'bed_file' must be provided")
  }
  
  # Parse input coordinates
  if (!is.null(bed_file)) {
    if (verbose) message("Reading loci from BED file: ", bed_file)
    loci_coords <- .parse_bed_file(bed_file, verbose)
  } else {
    if (verbose) message("Parsing ", length(loci), " loci coordinates")
    loci_coords <- .parse_loci_input(loci, verbose)
  }
  
  if (nrow(loci_coords) == 0) {
    if (verbose) message("No valid loci found")
    return(.create_empty_annotation_dataframe())
  }
  
  if (verbose) message("Retrieving annotations for ", nrow(loci_coords), " loci")
  
  # Create a temporary candidate dataset for these specific loci
  temp_candidates <- data.frame(
    chromosome = loci_coords$chromosome,
    position = loci_coords$position,
    stringsAsFactors = FALSE
  )
  
  # Use compile_funseq_results infrastructure to get comprehensive data
  # First try with enrichment data if requested
  if (include_enrichment) {
    tryCatch({
      # Get all available annotation data using compile_funseq_results
      all_annotations <- compile_funseq_results(
        con = con,
        stage = "annotations", 
        include = c("GO", "KEGG", "Pfam", "InterPro", "eggNOG"),
        candidate_loci = temp_candidates,
        verbose = FALSE
      )
      
      # Try to add enrichment information if available
      tryCatch({
        # Get available ORA analysis IDs
        ora_analyses <- DBI::dbGetQuery(con, "SELECT DISTINCT analysis_id FROM ora_results")
        if (nrow(ora_analyses) > 0) {
          result_data <- compile_funseq_results(
            con = con,
            stage = "enrichment",
            data = all_annotations,
            analysis_ids = ora_analyses$analysis_id,
            verbose = FALSE
          )
        } else {
          result_data <- all_annotations
          # Add empty enrichment columns
          result_data$enriched <- FALSE
          result_data$enrichment_fdr <- NA_real_
          result_data$enrichment_pvalue <- NA_real_
          result_data$enriched_terms <- ""
          result_data$enriched_term_ids <- ""
          result_data$enrichment_analysis_ids <- ""
        }
      }, error = function(e) {
        if (verbose) message("    - Enrichment data not available: ", e$message)
        result_data <<- all_annotations
        # Add empty enrichment columns
        result_data$enriched <<- FALSE
        result_data$enrichment_fdr <<- NA_real_
        result_data$enrichment_pvalue <<- NA_real_
        result_data$enriched_terms <<- ""
        result_data$enriched_term_ids <<- ""
        result_data$enrichment_analysis_ids <<- ""
      })
      
      # Filter results to only include the requested loci
      if (nrow(result_data) > 0) {
        result_data <- .filter_to_requested_loci(result_data, loci_coords, verbose)
      }
      
    }, error = function(e) {
      if (verbose) message("    - Error retrieving annotation data: ", e$message)
      result_data <<- .create_empty_annotation_dataframe()
    })
  } else {
    # Just get annotation data without enrichment
    tryCatch({
      result_data <- compile_funseq_results(
        con = con,
        stage = "annotations", 
        include = c("GO", "KEGG", "Pfam", "InterPro", "eggNOG"),
        candidate_loci = temp_candidates,
        verbose = FALSE
      )
      # Add empty enrichment columns
      result_data$enriched <- FALSE
      result_data$enrichment_fdr <- NA_real_
      result_data$enrichment_pvalue <- NA_real_
      result_data$enriched_terms <- ""
      result_data$enriched_term_ids <- ""
      result_data$enrichment_analysis_ids <- ""
      
      # Filter results to only include the requested loci
      if (nrow(result_data) > 0) {
        result_data <- .filter_to_requested_loci(result_data, loci_coords, verbose)
      }
    }, error = function(e) {
      if (verbose) message("    - Error retrieving annotation data: ", e$message)
      result_data <- .create_empty_annotation_dataframe()
    })
  }
  
  # Ensure all requested loci are represented, even if no annotations found
  result_data <- .ensure_all_loci_present(result_data, loci_coords, verbose)
  
  if (verbose) {
    annotated_count <- sum(!is.na(result_data$gene_name) & result_data$gene_name != "")
    enriched_count <- sum(result_data$enriched, na.rm = TRUE)
    message("Retrieved annotations for ", nrow(result_data), " loci:")
    message("  - ", annotated_count, " have functional annotations")
    if (include_enrichment) {
      message("  - ", enriched_count, " are associated with enriched terms")
    }
  }
  
  return(result_data)
}

# Helper functions for get_locus_annotations()

#' Parse loci input into standardized coordinate format
#' @keywords internal
.parse_loci_input <- function(loci, verbose) {
  if (is.data.frame(loci)) {
    # Input is already a data.frame
    if (all(c("chromosome", "position") %in% names(loci))) {
      result <- data.frame(
        chromosome = as.character(loci$chromosome),
        position = as.numeric(loci$position),
        stringsAsFactors = FALSE
      )
      # Remove rows with missing coordinates
      result <- result[!is.na(result$position), ]
      if (verbose) message("  - Parsed ", nrow(result), " loci from data.frame")
      return(result)
    } else {
      stop("Data.frame input must contain 'chromosome' and 'position' columns")
    }
  } else if (is.character(loci)) {
    # Parse character vector of coordinates
    result <- data.frame(
      chromosome = character(0),
      position = numeric(0),
      stringsAsFactors = FALSE
    )
    
    for (locus in loci) {
      # Handle different formats: "chr:pos", "chr pos", "chr_pos"
      if (grepl(":", locus)) {
        parts <- strsplit(locus, ":")[[1]]
      } else if (grepl("\\s+", locus)) {
        parts <- strsplit(locus, "\\s+")[[1]]
      } else if (grepl("_", locus)) {
        parts <- strsplit(locus, "_")[[1]]
      } else {
        if (verbose) message("  - Skipping invalid locus format: ", locus)
        next
      }
      
      if (length(parts) == 2) {
        chr <- trimws(parts[1])
        pos <- suppressWarnings(as.numeric(gsub("[^0-9]", "", parts[2])))
        
        if (!is.na(pos)) {
          result <- rbind(result, data.frame(
            chromosome = chr,
            position = pos,
            stringsAsFactors = FALSE
          ))
        } else {
          if (verbose) message("  - Skipping invalid position: ", locus)
        }
      } else {
        if (verbose) message("  - Skipping invalid format: ", locus)
      }
    }
    
    if (verbose) message("  - Parsed ", nrow(result), "/", length(loci), " valid loci")
    return(result)
  } else {
    stop("loci must be either a character vector or a data.frame")
  }
}

#' Parse BED file into coordinate format
#' @keywords internal  
.parse_bed_file <- function(bed_file, verbose) {
  if (!file.exists(bed_file)) {
    stop("BED file not found: ", bed_file)
  }
  
  # Read BED file (tab-separated, first 3 columns are chr, start, end)
  bed_data <- read.table(bed_file, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
  
  if (ncol(bed_data) < 3) {
    stop("BED file must have at least 3 columns (chromosome, start, end)")
  }
  
  # Use middle position of each interval
  result <- data.frame(
    chromosome = as.character(bed_data[, 1]),
    position = as.numeric(bed_data[, 2]) + floor((as.numeric(bed_data[, 3]) - as.numeric(bed_data[, 2])) / 2),
    stringsAsFactors = FALSE
  )
  
  # Remove rows with missing coordinates
  result <- result[!is.na(result$position), ]
  
  if (verbose) message("  - Parsed ", nrow(result), " regions from BED file")
  return(result)
}

#' Create empty annotation dataframe with all expected columns
#' @keywords internal
.create_empty_annotation_dataframe <- function() {
  data.frame(
    locus_id = character(0),
    chromosome = character(0),
    position = numeric(0),
    gene_name = character(0),
    protein_name = character(0),
    uniprot_accession = character(0),
    go_terms = character(0),
    go_names = character(0),
    go_categories = character(0),
    kegg_pathways = character(0),
    kegg_pathway_names = character(0),
    kegg_brite_categories = character(0),
    kegg_modules = character(0),
    pfam_domains = character(0),
    pfam_domain_names = character(0),
    interpro_families = character(0),
    interpro_family_names = character(0),
    eggnog_categories = character(0),
    eggnog_taxonomic_scopes = character(0),
    dataset_type = character(0),
    enriched = logical(0),
    enrichment_fdr = numeric(0),
    enrichment_pvalue = numeric(0),
    enriched_terms = character(0),
    enriched_term_ids = character(0),
    enrichment_analysis_ids = character(0),
    stringsAsFactors = FALSE
  )
}

#' Ensure all requested loci are present in results
#' @keywords internal
.ensure_all_loci_present <- function(result_data, loci_coords, verbose) {
  # Check which loci are missing from results
  missing_loci <- data.frame(
    chromosome = character(0),
    position = numeric(0),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:nrow(loci_coords)) {
    chr <- loci_coords$chromosome[i]
    pos <- loci_coords$position[i]
    
    # Check if this locus exists in results
    exists <- any(result_data$chromosome == chr & result_data$position == pos)
    
    if (!exists) {
      missing_loci <- rbind(missing_loci, data.frame(
        chromosome = chr,
        position = pos,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Add missing loci with NA values
  if (nrow(missing_loci) > 0) {
    if (verbose) message("  - Adding ", nrow(missing_loci), " loci with no annotations")
    
    # Create empty rows for missing loci
    empty_rows <- data.frame(
      locus_id = paste0(missing_loci$chromosome, ":", missing_loci$position),
      chromosome = missing_loci$chromosome,
      position = missing_loci$position,
      gene_name = NA_character_,
      protein_name = NA_character_,
      uniprot_accession = NA_character_,
      go_terms = NA_character_,
      go_names = NA_character_,
      go_categories = NA_character_,
      kegg_pathways = NA_character_,
      kegg_pathway_names = NA_character_,
      kegg_brite_categories = NA_character_,
      kegg_modules = NA_character_,
      pfam_domains = NA_character_,
      pfam_domain_names = NA_character_,
      interpro_families = NA_character_,
      interpro_family_names = NA_character_,
      eggnog_categories = NA_character_,
      eggnog_taxonomic_scopes = NA_character_,
      dataset_type = "candidate",
      enriched = FALSE,
      enrichment_fdr = NA_real_,
      enrichment_pvalue = NA_real_,
      enriched_terms = "",
      enriched_term_ids = "",
      enrichment_analysis_ids = "",
      stringsAsFactors = FALSE
    )
    
    # Add any missing columns that might exist in result_data
    for (col in names(result_data)) {
      if (!col %in% names(empty_rows)) {
        empty_rows[[col]] <- NA
      }
    }
    
    # Ensure column order matches
    empty_rows <- empty_rows[, names(result_data)]
    
    # Combine results
    result_data <- rbind(result_data, empty_rows)
  }
  
  # Sort by chromosome and position
  result_data <- result_data[order(result_data$chromosome, result_data$position), ]
  
  return(result_data)
}

#' Filter comprehensive results to only include requested loci
#' @keywords internal
.filter_to_requested_loci <- function(result_data, loci_coords, verbose) {
  # Create a matching key for requested loci
  requested_keys <- paste0(loci_coords$chromosome, ":", loci_coords$position)
  
  # Create matching key for result data
  result_keys <- paste0(result_data$chromosome, ":", result_data$position)
  
  # Filter to only include requested loci
  filtered_data <- result_data[result_keys %in% requested_keys, ]
  
  if (verbose && nrow(filtered_data) < nrow(result_data)) {
    message("    - Filtered from ", nrow(result_data), " to ", nrow(filtered_data), " loci")
  }
  
  return(filtered_data)
}


