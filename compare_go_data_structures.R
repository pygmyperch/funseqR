# Compare GO Data Structures Between Old and New Functions

compare_go_data_structures <- function(annotations, candidate_loci) {
  
  cat("=== Comparing GO Data Structures ===\n\n")
  
  # 1. Get candidate loci using new function logic
  has_chrom <- "chrom" %in% colnames(candidate_loci)
  candidate_df <- data.frame(
    chromosome = candidate_loci[[if (has_chrom) "chrom" else "chromosome"]],
    position = candidate_loci$start + 1,
    stringsAsFactors = FALSE
  )
  
  candidate_patterns <- with(candidate_df, paste0("_", chromosome, "_", position, "$"))
  candidate_locus_ids <- annotations$locus_id[
    sapply(annotations$locus_id, function(locus_id) {
      any(sapply(candidate_patterns, function(pattern) {
        grepl(pattern, locus_id)
      }))
    })
  ]
  
  cat("Candidate loci found:", length(candidate_locus_ids), "\n\n")
  
  # 2. Replicate NEW function GO data preparation
  cat("=== NEW FUNCTION APPROACH ===\n")
  
  go_annotations <- annotations[!is.na(annotations$go_terms) & annotations$go_terms != "", ]
  go_term_map_new <- list()
  
  for (i in 1:nrow(go_annotations)) {
    locus_id <- go_annotations$locus_id[i]
    go_terms <- unlist(strsplit(go_annotations$go_terms[i], ";"))
    go_categories <- unlist(strsplit(go_annotations$go_categories[i], ";"))
    
    for (j in seq_along(go_terms)) {
      go_id <- go_terms[j]
      go_category <- if (j <= length(go_categories)) go_categories[j] else "Unknown"
      
      if (!is.na(go_id) && go_id != "") {
        if (!go_id %in% names(go_term_map_new)) {
          go_term_map_new[[go_id]] <- list(
            go_id = go_id,
            go_category = go_category,
            loci = character(0)
          )
        }
        go_term_map_new[[go_id]]$loci <- c(go_term_map_new[[go_id]]$loci, locus_id)
      }
    }
  }
  
  # Filter for BP terms
  bp_terms_new <- go_term_map_new[sapply(go_term_map_new, function(x) x$go_category == "P")]
  bp_terms_sized_new <- bp_terms_new[sapply(bp_terms_new, function(x) length(unique(x$loci)) >= 3)]
  
  cat("Total GO terms (new):", length(go_term_map_new), "\n")
  cat("BP terms (new):", length(bp_terms_new), "\n")
  cat("BP terms sized >=3 (new):", length(bp_terms_sized_new), "\n")
  
  # Count terms with candidates
  bp_terms_with_candidates_new <- 0
  for (term_info in bp_terms_sized_new) {
    if (length(intersect(candidate_locus_ids, term_info$loci)) > 0) {
      bp_terms_with_candidates_new <- bp_terms_with_candidates_new + 1
    }
  }
  cat("BP terms with candidates (new):", bp_terms_with_candidates_new, "\n\n")
  
  # 3. Simulate OLD function approach (manual test)
  cat("=== OLD FUNCTION APPROACH (Manual Test Logic) ===\n")
  
  # The manual test used a different approach - let's replicate it
  go_term_map_old <- list()
  
  # Use only first 100 annotations like manual test did
  for (i in 1:min(nrow(go_annotations), 100)) {
    locus_id <- go_annotations$locus_id[i]
    
    if (!is.na(go_annotations$go_terms[i]) && go_annotations$go_terms[i] != "") {
      go_terms <- unlist(strsplit(go_annotations$go_terms[i], ";"))
      go_categories <- unlist(strsplit(go_annotations$go_categories[i], ";"))
      
      for (j in seq_along(go_terms)) {
        go_id <- go_terms[j]
        go_category <- if (j <= length(go_categories)) go_categories[j] else "Unknown"
        
        if (!is.na(go_id) && go_id != "") {
          if (!go_id %in% names(go_term_map_old)) {
            go_term_map_old[[go_id]] <- list(
              go_id = go_id,
              go_category = go_category,
              loci = character(0)
            )
          }
          go_term_map_old[[go_id]]$loci <- c(go_term_map_old[[go_id]]$loci, locus_id)
        }
      }
    }
  }
  
  # Filter for BP terms
  bp_terms_old <- go_term_map_old[sapply(go_term_map_old, function(x) x$go_category == "P")]
  bp_terms_sized_old <- bp_terms_old[sapply(bp_terms_old, function(x) length(unique(x$loci)) >= 3)]
  
  cat("Total GO terms (old/limited):", length(go_term_map_old), "\n")
  cat("BP terms (old/limited):", length(bp_terms_old), "\n")
  cat("BP terms sized >=3 (old/limited):", length(bp_terms_sized_old), "\n")
  
  # Count terms with candidates
  bp_terms_with_candidates_old <- 0
  for (term_info in bp_terms_sized_old) {
    if (length(intersect(candidate_locus_ids, term_info$loci)) > 0) {
      bp_terms_with_candidates_old <- bp_terms_with_candidates_old + 1
    }
  }
  cat("BP terms with candidates (old/limited):", bp_terms_with_candidates_old, "\n\n")
  
  # 4. Key Differences Analysis
  cat("=== KEY DIFFERENCES ===\n")
  cat("The manual test limitation explains the discrepancy:\n")
  cat("  - Manual test: processed only 100 annotations\n")
  cat("  - New function: processes all", nrow(go_annotations), "annotations\n")
  cat("  - This creates different GO term sets and candidate overlaps\n\n")
  
  # 5. Show some example terms that differ
  cat("=== SAMPLE TERM ANALYSIS ===\n")
  
  # Get a few terms from each and compare
  sample_new <- names(bp_terms_sized_new)[1:min(5, length(bp_terms_sized_new))]
  sample_old <- names(bp_terms_sized_old)[1:min(5, length(bp_terms_sized_old))]
  
  cat("Sample NEW terms:", paste(sample_new, collapse = ", "), "\n")
  cat("Sample OLD terms:", paste(sample_old, collapse = ", "), "\n")
  
  # Check overlap between the two sets
  common_terms <- intersect(sample_new, sample_old)
  cat("Common terms:", paste(common_terms, collapse = ", "), "\n")
  
  return(list(
    new_bp_count = length(bp_terms_sized_new),
    old_bp_count = length(bp_terms_sized_old),
    new_with_candidates = bp_terms_with_candidates_new,
    old_with_candidates = bp_terms_with_candidates_old
  ))
}

# Usage:
# comparison <- compare_go_data_structures(annotations, candidates$bed_file)