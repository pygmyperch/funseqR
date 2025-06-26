# Debug Term Matching in New vs Old Functions
# This script investigates why candidate-term matching differs between functions

debug_term_matching <- function(annotations, candidate_loci, ontology = "BP") {
  
  cat("=== Debugging Term Matching ===\n")
  cat("Ontology:", ontology, "\n\n")
  
  # Step 1: Get candidate loci (replicate new function logic)
  cat("Step 1: Candidate Loci Processing\n")
  
  has_chrom <- "chrom" %in% colnames(candidate_loci)
  has_chromosome <- "chromosome" %in% colnames(candidate_loci)
  has_start <- "start" %in% colnames(candidate_loci)
  has_end <- "end" %in% colnames(candidate_loci)
  
  if ((has_chrom || has_chromosome) && has_start && has_end) {
    chrom_col <- if (has_chrom) "chrom" else "chromosome"
    candidate_df <- data.frame(
      chromosome = candidate_loci[[chrom_col]],
      position = candidate_loci$start + 1,  # Convert 0-based to 1-based
      stringsAsFactors = FALSE
    )
  }
  
  candidate_patterns <- with(candidate_df, paste0("_", chromosome, "_", position, "$"))
  candidate_locus_ids <- annotations$locus_id[
    sapply(annotations$locus_id, function(locus_id) {
      any(sapply(candidate_patterns, function(pattern) {
        grepl(pattern, locus_id)
      }))
    })
  ]
  
  cat("  Candidate locus IDs found:", length(candidate_locus_ids), "\n")
  cat("  First few candidate IDs:", paste(head(candidate_locus_ids, 3), collapse = ", "), "\n\n")
  
  # Step 2: Process GO data (replicate new function logic)
  cat("Step 2: GO Data Processing (New Function Logic)\n")
  
  go_annotations <- annotations[!is.na(annotations$go_terms) & annotations$go_terms != "", ]
  cat("  GO annotations count:", nrow(go_annotations), "\n")
  
  # Build GO term map like new function
  go_term_map <- list()
  
  for (i in 1:nrow(go_annotations)) {
    locus_id <- go_annotations$locus_id[i]
    go_terms <- unlist(strsplit(go_annotations$go_terms[i], ";"))
    go_categories <- unlist(strsplit(go_annotations$go_categories[i], ";"))
    
    for (j in seq_along(go_terms)) {
      go_id <- go_terms[j]
      go_category <- if (j <= length(go_categories)) go_categories[j] else "Unknown"
      
      if (!is.na(go_id) && go_id != "") {
        if (!go_id %in% names(go_term_map)) {
          go_term_map[[go_id]] <- list(
            go_id = go_id,
            go_category = go_category,
            loci = character(0)
          )
        }
        go_term_map[[go_id]]$loci <- c(go_term_map[[go_id]]$loci, locus_id)
      }
    }
  }
  
  cat("  Total GO terms built:", length(go_term_map), "\n")
  
  # Filter for BP terms
  ontology_map <- c("BP" = "P", "MF" = "F", "CC" = "C")
  target_category <- ontology_map[ontology]
  
  bp_terms <- go_term_map[sapply(go_term_map, function(x) x$go_category == target_category)]
  cat("  BP terms available:", length(bp_terms), "\n")
  
  # Filter by size (min_genes = 3)
  bp_terms_sized <- bp_terms[sapply(bp_terms, function(x) length(unique(x$loci)) >= 3)]
  cat("  BP terms passing size filter (>=3):", length(bp_terms_sized), "\n")
  
  # Check candidate overlap for each term
  cat("\nStep 3: Candidate-Term Overlap Analysis\n")
  
  terms_with_candidates <- 0
  terms_without_candidates <- 0
  
  overlap_details <- list()
  
  for (i in 1:min(10, length(bp_terms_sized))) {  # Check first 10 terms
    term_info <- bp_terms_sized[[i]]
    go_id <- term_info$go_id
    term_loci <- unique(term_info$loci)
    
    # Check overlap with candidates
    candidate_overlap <- intersect(candidate_locus_ids, term_loci)
    overlap_count <- length(candidate_overlap)
    
    overlap_details[[go_id]] <- list(
      go_id = go_id,
      total_loci = length(term_loci),
      candidate_overlap = overlap_count,
      candidate_overlap_ids = if(overlap_count > 0) candidate_overlap else NULL
    )
    
    if (overlap_count > 0) {
      terms_with_candidates <- terms_with_candidates + 1
      cat("  ", go_id, ": ", overlap_count, " candidates /", length(term_loci), " background\n")
      if (overlap_count <= 3) {
        cat("    Candidate IDs:", paste(candidate_overlap, collapse = ", "), "\n")
      }
    } else {
      terms_without_candidates <- terms_without_candidates + 1
    }
  }
  
  cat("\nStep 4: Summary\n")
  cat("  Terms with candidates (first 10):", terms_with_candidates, "\n")
  cat("  Terms without candidates (first 10):", terms_without_candidates, "\n")
  
  # Step 5: Compare with old function approach
  cat("\nStep 5: Comparison Check\n")
  cat("  Manual test found 6 testable BP terms with candidates\n")
  cat("  New function reports 13 BP terms tested\n")
  cat("  Discrepancy suggests different candidate-term matching logic\n")
  
  return(list(
    candidate_locus_ids = candidate_locus_ids,
    bp_terms_count = length(bp_terms),
    bp_terms_sized_count = length(bp_terms_sized),
    overlap_details = overlap_details
  ))
}

# Usage:
# debug_results <- debug_term_matching(annotations, candidates$bed_file, "BP")