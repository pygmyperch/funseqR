#' Run GO Enrichment Analysis on Processed Annotations
#'
#' Performs GO enrichment analysis using candidate loci and processed annotation data.
#' This is a simplified interface that works with the standardized output from process_annotations().
#'
#' @param annotations Data frame from process_annotations() containing GO annotations
#' @param candidate_loci Data frame with coordinates (BED or VCF format), or vector of locus_ids. 
#'   Accepts: (1) BED format with 'chrom' (or 'chromosome'), 'start', 'end' columns, 
#'   (2) VCF format with 'chrom' (or 'chromosome'), 'position' columns, or (3) vector of locus_ids
#' @param ontologies Character vector. GO ontologies to test: c("BP", "MF", "CC"). Default is c("BP", "MF", "CC")
#' @param min_genes Integer. Minimum genes for GO term testing. Default is 5
#' @param max_genes Integer. Maximum genes for GO term testing. Default is 500
#' @param significance_threshold Numeric. FDR threshold for significance. Default is 0.05
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return List containing enrichment results for each ontology tested
#'
#' @examples
#' \dontrun{
#' # Process annotations first
#' annotations <- process_annotations(con, include = c("GO", "KEGG"))
#' 
#' # Using VCF format
#' candidates_vcf <- data.frame(chromosome = c("LG1", "LG2"), position = c(12345, 67890))
#' go_results <- run_go_enrichment_analysis(annotations, candidates_vcf)
#' 
#' # Using BED format (from import_candidate_loci)
#' candidates <- import_candidate_loci(con, candidate_vcf_file, background_file_id = 1)
#' go_results <- run_go_enrichment_analysis(annotations, candidates$bed_file)
#' 
#' # Using locus IDs directly
#' go_results <- run_go_enrichment_analysis(annotations, c("1_LG1_12345", "1_LG2_67890"))
#' 
#' print(go_results$BP)
#' }
#'
#' @export
run_go_enrichment_analysis <- function(annotations, candidate_loci, 
                                     ontologies = c("BP", "MF", "CC"),
                                     min_genes = 5, max_genes = 500,
                                     significance_threshold = 0.05,
                                     verbose = TRUE) {
  
  if (verbose) message("Running GO enrichment analysis...")
  
  # Check if GO annotations are available
  if (!"go_terms" %in% colnames(annotations)) {
    stop("GO annotations not found. Run process_annotations() with include='GO' first.")
  }
  
  # Identify candidate loci - handle multiple input formats with flexible column naming
  if (is.data.frame(candidate_loci)) {
    
    # Detect available column names (support both standard formats)
    has_chrom <- "chrom" %in% colnames(candidate_loci)
    has_chromosome <- "chromosome" %in% colnames(candidate_loci)
    has_start <- "start" %in% colnames(candidate_loci)
    has_end <- "end" %in% colnames(candidate_loci)
    has_position <- "position" %in% colnames(candidate_loci)
    
    # Handle BED file format (chrom/chromosome, start, end)
    if ((has_chrom || has_chromosome) && has_start && has_end) {
      chrom_col <- if (has_chrom) "chrom" else "chromosome"
      if (verbose) message("  - Detected BED format input (", chrom_col, ", start, end)")
      
      # Convert BED coordinates - add 1 to start position to convert from 0-based to 1-based
      candidate_df <- data.frame(
        chromosome = candidate_loci[[chrom_col]],
        position = candidate_loci$start + 1,  # Convert 0-based BED to 1-based VCF coordinates
        stringsAsFactors = FALSE
      )
    } 
    # Handle VCF format (chrom/chromosome, position)
    else if ((has_chrom || has_chromosome) && has_position) {
      chrom_col <- if (has_chrom) "chrom" else "chromosome"
      if (verbose) message("  - Detected VCF format input (", chrom_col, ", position)")
      
      candidate_df <- data.frame(
        chromosome = candidate_loci[[chrom_col]],
        position = candidate_loci$position,
        stringsAsFactors = FALSE
      )
    }
    # Handle file_id input - extract from database
    else if ("file_id" %in% colnames(candidate_loci) && nrow(candidate_loci) == 1) {
      if (verbose) message("  - Extracting candidate loci from file_id: ", candidate_loci$file_id[1])
      candidate_df <- DBI::dbGetQuery(annotations_con, 
        "SELECT chromosome, position FROM vcf_data WHERE file_id = ?",
        list(candidate_loci$file_id[1]))
    }
    else {
      # Show helpful error with actual column names found
      found_cols <- paste(colnames(candidate_loci), collapse = ", ")
      stop("candidate_loci data frame format not recognized.\n",
           "Found columns: ", found_cols, "\n",
           "Expected one of:\n",
           "  - BED format: 'chrom' (or 'chromosome'), 'start', 'end'\n",
           "  - VCF format: 'chrom' (or 'chromosome'), 'position'\n",
           "  - File reference: 'file_id'")
    }
    
    # Create coordinate patterns for matching
    # Annotation locus_ids have format: vcf_id_chromosome_position
    # We need to match on chromosome_position part
    candidate_patterns <- with(candidate_df, paste0("_", chromosome, "_", position, "$"))
    
    # Match against locus_id patterns in annotations
    candidate_locus_ids <- annotations$locus_id[
      sapply(annotations$locus_id, function(locus_id) {
        any(sapply(candidate_patterns, function(pattern) {
          grepl(pattern, locus_id)
        }))
      })
    ]
  } 
  # Handle direct file_id input
  else if (is.numeric(candidate_loci) && length(candidate_loci) == 1) {
    if (verbose) message("  - Extracting candidate loci from file_id: ", candidate_loci)
    # Need to access the database connection - this requires the con parameter
    stop("Direct file_id input not yet supported. Please provide a data frame with coordinates.")
  }
  else {
    # Assume vector of locus_ids
    candidate_locus_ids <- candidate_loci
  }
  
  if (length(candidate_locus_ids) == 0) {
    if (verbose && exists("candidate_patterns")) {
      message("  - Debug: First few candidate patterns: ", paste(head(candidate_patterns, 3), collapse = ", "))
      message("  - Debug: First few annotation locus_ids: ", paste(head(annotations$locus_id, 3), collapse = ", "))
    }
    stop("No candidate loci found in annotations")
  }
  
  if (verbose) {
    message("  - Candidate loci: ", length(candidate_locus_ids))
    message("  - Background loci: ", nrow(annotations))
    message("  - Testing ontologies: ", paste(ontologies, collapse = ", "))
  }
  
  # Filter annotations to only those with GO terms
  go_annotations <- annotations[!is.na(annotations$go_terms) & annotations$go_terms != "", ]
  
  if (nrow(go_annotations) == 0) {
    stop("No GO annotations found in the data")
  }
  
  # Prepare GO data
  go_data <- .prepare_go_data_from_annotations(go_annotations, candidate_locus_ids, verbose)
  
  # Run enrichment for each ontology
  results <- list()
  
  for (ontology in ontologies) {
    if (verbose) message("  - Analyzing ", ontology, " ontology...")
    
    ontology_results <- .perform_go_enrichment_from_data(
      go_data, ontology, min_genes, max_genes, significance_threshold, verbose
    )
    
    results[[ontology]] <- ontology_results
  }
  
  if (verbose) {
    total_significant <- sum(sapply(results, function(x) sum(x$p_adjusted < significance_threshold, na.rm = TRUE)))
    message("  - Total significant terms found: ", total_significant)
  }
  
  return(results)
}

#' Run KEGG Pathway Enrichment Analysis on Processed Annotations
#'
#' Performs KEGG pathway enrichment analysis using candidate loci and processed annotation data.
#'
#' @param annotations Data frame from process_annotations() containing KEGG annotations
#' @param candidate_loci Data frame with coordinates (BED or VCF format), or vector of locus_ids. 
#'   Accepts: (1) BED format with 'chrom' (or 'chromosome'), 'start', 'end' columns, 
#'   (2) VCF format with 'chrom' (or 'chromosome'), 'position' columns, or (3) vector of locus_ids
#' @param min_pathways Integer. Minimum genes for pathway testing. Default is 3
#' @param max_pathways Integer. Maximum genes for pathway testing. Default is 500
#' @param significance_threshold Numeric. FDR threshold for significance. Default is 0.05
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return Data frame with KEGG pathway enrichment results
#'
#' @examples
#' \dontrun{
#' # Process annotations first
#' annotations <- process_annotations(con, include = c("GO", "KEGG"))
#' 
#' # Using BED format (from import_candidate_loci) 
#' candidates <- import_candidate_loci(con, candidate_vcf_file, background_file_id = 1)
#' kegg_results <- run_kegg_enrichment_analysis(annotations, candidates$bed_file)
#' 
#' # Using VCF format
#' candidates_vcf <- data.frame(chromosome = c("LG1", "LG2"), position = c(12345, 67890))
#' kegg_results <- run_kegg_enrichment_analysis(annotations, candidates_vcf)
#' 
#' print(kegg_results)
#' }
#'
#' @export
run_kegg_enrichment_analysis <- function(annotations, candidate_loci,
                                       min_pathways = 3, max_pathways = 500,
                                       significance_threshold = 0.05,
                                       verbose = TRUE) {
  
  if (verbose) message("Running KEGG pathway enrichment analysis...")
  
  # Check if KEGG annotations are available
  if (!"kegg_pathways" %in% colnames(annotations)) {
    stop("KEGG annotations not found. Run process_annotations() with include='KEGG' first.")
  }
  
  # Identify candidate loci - handle multiple input formats with flexible column naming
  if (is.data.frame(candidate_loci)) {
    
    # Detect available column names (support both standard formats)
    has_chrom <- "chrom" %in% colnames(candidate_loci)
    has_chromosome <- "chromosome" %in% colnames(candidate_loci)
    has_start <- "start" %in% colnames(candidate_loci)
    has_end <- "end" %in% colnames(candidate_loci)
    has_position <- "position" %in% colnames(candidate_loci)
    
    # Handle BED file format (chrom/chromosome, start, end)
    if ((has_chrom || has_chromosome) && has_start && has_end) {
      chrom_col <- if (has_chrom) "chrom" else "chromosome"
      if (verbose) message("  - Detected BED format input (", chrom_col, ", start, end)")
      
      # Convert BED coordinates - add 1 to start position to convert from 0-based to 1-based
      candidate_df <- data.frame(
        chromosome = candidate_loci[[chrom_col]],
        position = candidate_loci$start + 1,  # Convert 0-based BED to 1-based VCF coordinates
        stringsAsFactors = FALSE
      )
    } 
    # Handle VCF format (chrom/chromosome, position)
    else if ((has_chrom || has_chromosome) && has_position) {
      chrom_col <- if (has_chrom) "chrom" else "chromosome"
      if (verbose) message("  - Detected VCF format input (", chrom_col, ", position)")
      
      candidate_df <- data.frame(
        chromosome = candidate_loci[[chrom_col]],
        position = candidate_loci$position,
        stringsAsFactors = FALSE
      )
    }
    else {
      # Show helpful error with actual column names found
      found_cols <- paste(colnames(candidate_loci), collapse = ", ")
      stop("candidate_loci data frame format not recognized.\n",
           "Found columns: ", found_cols, "\n",
           "Expected one of:\n",
           "  - BED format: 'chrom' (or 'chromosome'), 'start', 'end'\n",
           "  - VCF format: 'chrom' (or 'chromosome'), 'position'\n",
           "  - File reference: 'file_id'")
    }
    
    # Create coordinate patterns for matching
    # Annotation locus_ids have format: vcf_id_chromosome_position
    # We need to match on chromosome_position part
    candidate_patterns <- with(candidate_df, paste0("_", chromosome, "_", position, "$"))
    
    # Match against locus_id patterns in annotations
    candidate_locus_ids <- annotations$locus_id[
      sapply(annotations$locus_id, function(locus_id) {
        any(sapply(candidate_patterns, function(pattern) {
          grepl(pattern, locus_id)
        }))
      })
    ]
  } else {
    candidate_locus_ids <- candidate_loci
  }
  
  if (length(candidate_locus_ids) == 0) {
    if (verbose && exists("candidate_patterns")) {
      message("  - Debug: First few candidate patterns: ", paste(head(candidate_patterns, 3), collapse = ", "))
      message("  - Debug: First few annotation locus_ids: ", paste(head(annotations$locus_id, 3), collapse = ", "))
    }
    stop("No candidate loci found in annotations")
  }
  
  # Filter annotations to only those with KEGG pathways
  kegg_annotations <- annotations[!is.na(annotations$kegg_pathways) & annotations$kegg_pathways != "", ]
  
  if (nrow(kegg_annotations) == 0) {
    stop("No KEGG annotations found in the data")
  }
  
  if (verbose) {
    message("  - Candidate loci: ", length(candidate_locus_ids))
    message("  - Background loci with KEGG: ", nrow(kegg_annotations))
  }
  
  # Prepare KEGG pathway data
  kegg_data <- .prepare_kegg_data_from_annotations(kegg_annotations, candidate_locus_ids, verbose)
  
  # Perform enrichment analysis
  results <- .perform_kegg_enrichment_from_data(
    kegg_data, min_pathways, max_pathways, significance_threshold, verbose
  )
  
  if (verbose) {
    significant_count <- sum(results$p_adjusted < significance_threshold, na.rm = TRUE)
    message("  - Significant pathways found: ", significant_count)
  }
  
  return(results)
}

#' Prepare GO data from processed annotations
#' @keywords internal
.prepare_go_data_from_annotations <- function(annotations, candidate_locus_ids, verbose) {
  
  # Create GO term mappings
  go_term_map <- list()
  
  for (i in 1:nrow(annotations)) {
    locus_id <- annotations$locus_id[i]
    go_terms <- unlist(strsplit(annotations$go_terms[i], ";"))
    go_names <- unlist(strsplit(annotations$go_names[i], ";"))
    go_categories <- unlist(strsplit(annotations$go_categories[i], ";"))
    
    for (j in seq_along(go_terms)) {
      go_id <- go_terms[j]
      if (!is.na(go_id) && go_id != "") {
        if (!go_id %in% names(go_term_map)) {
          go_term_map[[go_id]] <- list(
            go_id = go_id,
            go_name = if (j <= length(go_names)) go_names[j] else "Unknown",
            go_category = if (j <= length(go_categories)) go_categories[j] else "Unknown",
            loci = character(0)
          )
        }
        go_term_map[[go_id]]$loci <- c(go_term_map[[go_id]]$loci, locus_id)
      }
    }
  }
  
  # Convert to data frame format expected by enrichment functions
  go_data <- list(
    foreground = list(
      loci = candidate_locus_ids,
      genes = candidate_locus_ids  # Using locus_ids as gene identifiers
    ),
    background = list(
      loci = annotations$locus_id,
      genes = annotations$locus_id
    ),
    go_terms = do.call(rbind, lapply(go_term_map, function(x) {
      data.frame(
        go_id = x$go_id,
        go_name = x$go_name,
        go_category = x$go_category,
        locus_count = length(unique(x$loci)),
        stringsAsFactors = FALSE
      )
    })),
    term_mappings = go_term_map
  )
  
  if (verbose) {
    message("    - GO terms found: ", length(go_term_map))
    message("    - Candidate loci with GO: ", 
            length(intersect(candidate_locus_ids, annotations$locus_id)))
  }
  
  return(go_data)
}

#' Prepare KEGG data from processed annotations  
#' @keywords internal
.prepare_kegg_data_from_annotations <- function(annotations, candidate_locus_ids, verbose) {
  
  # Create KEGG pathway mappings
  pathway_map <- list()
  
  for (i in 1:nrow(annotations)) {
    locus_id <- annotations$locus_id[i]
    pathways <- unlist(strsplit(annotations$kegg_pathways[i], ";"))
    pathway_names <- unlist(strsplit(annotations$kegg_pathway_names[i], ";"))
    
    for (j in seq_along(pathways)) {
      pathway_id <- pathways[j]
      if (!is.na(pathway_id) && pathway_id != "") {
        if (!pathway_id %in% names(pathway_map)) {
          pathway_map[[pathway_id]] <- list(
            pathway_id = pathway_id,
            pathway_name = if (j <= length(pathway_names)) pathway_names[j] else "Unknown",
            loci = character(0)
          )
        }
        pathway_map[[pathway_id]]$loci <- c(pathway_map[[pathway_id]]$loci, locus_id)
      }
    }
  }
  
  # Convert to expected format
  kegg_data <- list(
    foreground = list(
      loci = candidate_locus_ids,
      genes = candidate_locus_ids
    ),
    background = list(
      loci = annotations$locus_id,
      genes = annotations$locus_id
    ),
    pathways = do.call(rbind, lapply(pathway_map, function(x) {
      data.frame(
        pathway_id = x$pathway_id,
        pathway_name = x$pathway_name,
        locus_count = length(unique(x$loci)),
        stringsAsFactors = FALSE
      )
    })),
    pathway_mappings = pathway_map
  )
  
  if (verbose) {
    message("    - KEGG pathways found: ", length(pathway_map))
    message("    - Candidate loci with KEGG: ", 
            length(intersect(candidate_locus_ids, annotations$locus_id)))
  }
  
  return(kegg_data)
}

#' Perform GO enrichment analysis from prepared data
#' @keywords internal
.perform_go_enrichment_from_data <- function(go_data, ontology, min_genes, max_genes, 
                                           significance_threshold, verbose) {
  
  # Filter GO terms for the specific ontology
  ontology_terms <- go_data$go_terms[go_data$go_terms$go_category == ontology, ]
  
  if (nrow(ontology_terms) == 0) {
    if (verbose) message("    - No ", ontology, " terms found")
    return(data.frame())
  }
  
  # Filter by size
  size_filtered <- ontology_terms[ontology_terms$locus_count >= min_genes & 
                                 ontology_terms$locus_count <= max_genes, ]
  
  if (nrow(size_filtered) == 0) {
    if (verbose) message("    - No ", ontology, " terms pass size filters")
    return(data.frame())
  }
  
  # Perform enrichment tests
  results <- data.frame()
  
  for (i in 1:nrow(size_filtered)) {
    go_id <- size_filtered$go_id[i]
    term_loci <- go_data$term_mappings[[go_id]]$loci
    
    # Fisher's exact test
    candidate_with_term <- length(intersect(go_data$foreground$loci, term_loci))
    candidate_without_term <- length(go_data$foreground$loci) - candidate_with_term
    background_with_term <- length(intersect(go_data$background$loci, term_loci)) - candidate_with_term
    background_without_term <- length(go_data$background$loci) - length(intersect(go_data$background$loci, term_loci)) - candidate_without_term
    
    # Create contingency table
    cont_table <- matrix(c(candidate_with_term, candidate_without_term,
                          background_with_term, background_without_term), 
                        nrow = 2, byrow = TRUE)
    
    if (candidate_with_term > 0) {
      test_result <- fisher.test(cont_table, alternative = "greater")
      
      results <- rbind(results, data.frame(
        go_id = go_id,
        go_name = size_filtered$go_name[i],
        go_category = ontology,
        candidate_count = candidate_with_term,
        background_count = length(intersect(go_data$background$loci, term_loci)),
        total_candidates = length(go_data$foreground$loci),
        total_background = length(go_data$background$loci),
        fold_enrichment = (candidate_with_term / length(go_data$foreground$loci)) / 
                         (length(intersect(go_data$background$loci, term_loci)) / length(go_data$background$loci)),
        p_value = test_result$p.value,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  if (nrow(results) > 0) {
    # Adjust p-values
    results$p_adjusted <- p.adjust(results$p_value, method = "fdr")
    
    # Sort by p-value
    results <- results[order(results$p_value), ]
    rownames(results) <- NULL
  }
  
  return(results)
}

#' Perform KEGG enrichment analysis from prepared data
#' @keywords internal
.perform_kegg_enrichment_from_data <- function(kegg_data, min_pathways, max_pathways, 
                                             significance_threshold, verbose) {
  
  # Filter pathways by size
  size_filtered <- kegg_data$pathways[kegg_data$pathways$locus_count >= min_pathways & 
                                     kegg_data$pathways$locus_count <= max_pathways, ]
  
  if (nrow(size_filtered) == 0) {
    if (verbose) message("    - No pathways pass size filters")
    return(data.frame())
  }
  
  # Perform enrichment tests
  results <- data.frame()
  
  for (i in 1:nrow(size_filtered)) {
    pathway_id <- size_filtered$pathway_id[i]
    pathway_loci <- kegg_data$pathway_mappings[[pathway_id]]$loci
    
    # Fisher's exact test
    candidate_with_pathway <- length(intersect(kegg_data$foreground$loci, pathway_loci))
    candidate_without_pathway <- length(kegg_data$foreground$loci) - candidate_with_pathway
    background_with_pathway <- length(intersect(kegg_data$background$loci, pathway_loci)) - candidate_with_pathway
    background_without_pathway <- length(kegg_data$background$loci) - length(intersect(kegg_data$background$loci, pathway_loci)) - candidate_without_pathway
    
    # Create contingency table
    cont_table <- matrix(c(candidate_with_pathway, candidate_without_pathway,
                          background_with_pathway, background_without_pathway), 
                        nrow = 2, byrow = TRUE)
    
    if (candidate_with_pathway > 0) {
      test_result <- fisher.test(cont_table, alternative = "greater")
      
      results <- rbind(results, data.frame(
        pathway_id = pathway_id,
        pathway_name = size_filtered$pathway_name[i],
        candidate_count = candidate_with_pathway,
        background_count = length(intersect(kegg_data$background$loci, pathway_loci)),
        total_candidates = length(kegg_data$foreground$loci),
        total_background = length(kegg_data$background$loci),
        fold_enrichment = (candidate_with_pathway / length(kegg_data$foreground$loci)) / 
                         (length(intersect(kegg_data$background$loci, pathway_loci)) / length(kegg_data$background$loci)),
        p_value = test_result$p.value,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  if (nrow(results) > 0) {
    # Adjust p-values
    results$p_adjusted <- p.adjust(results$p_value, method = "fdr")
    
    # Sort by p-value
    results <- results[order(results$p_value), ]
    rownames(results) <- NULL
  }
  
  return(results)
}