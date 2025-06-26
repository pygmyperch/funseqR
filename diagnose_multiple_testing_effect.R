# Diagnose Multiple Testing Effect
# This script demonstrates how testing more terms affects FDR adjustment

diagnose_multiple_testing_effect <- function(annotations, candidate_loci, ontology = "BP") {
  
  cat("=== Diagnosing Multiple Testing Effect ===\n")
  cat("Testing how the number of terms affects FDR adjustment\n\n")
  
  # Get candidate loci using same logic as new function
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
  
  cat("Candidate loci found:", length(candidate_locus_ids), "\n")
  
  # Build GO term map (same as new function)
  go_annotations <- annotations[!is.na(annotations$go_terms) & annotations$go_terms != "", ]
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
  
  # Filter for BP terms and size
  bp_terms <- go_term_map[sapply(go_term_map, function(x) x$go_category == "P")]
  bp_terms_sized <- bp_terms[sapply(bp_terms, function(x) length(unique(x$loci)) >= 3 & length(unique(x$loci)) <= 500)]
  
  cat("BP terms passing size filter:", length(bp_terms_sized), "\n")
  
  # Function to run Fisher's test on a subset of terms
  test_subset <- function(terms_to_test, subset_name) {
    cat("\n--- Testing", subset_name, "---\n")
    cat("Number of terms to test:", length(terms_to_test), "\n")
    
    results <- data.frame()
    
    for (i in 1:length(terms_to_test)) {
      term_info <- terms_to_test[[i]]
      go_id <- term_info$go_id
      term_loci <- unique(term_info$loci)
      
      # Fisher's exact test (same logic as new function)
      candidate_with_term <- length(intersect(candidate_locus_ids, term_loci))
      candidate_without_term <- length(candidate_locus_ids) - candidate_with_term
      
      total_background_with_term <- length(intersect(go_annotations$locus_id, term_loci))
      background_only_with_term <- total_background_with_term - candidate_with_term
      background_only_without_term <- length(go_annotations$locus_id) - total_background_with_term - candidate_without_term
      
      if (candidate_with_term > 0) {
        cont_table <- matrix(c(candidate_with_term, candidate_without_term,
                              background_only_with_term, background_only_without_term), 
                            nrow = 2, byrow = TRUE)
        
        test_result <- fisher.test(cont_table, alternative = "greater")
        
        results <- rbind(results, data.frame(
          go_id = go_id,
          candidate_count = candidate_with_term,
          background_count = total_background_with_term,
          total_candidates = length(candidate_locus_ids),
          total_background = length(go_annotations$locus_id),
          fold_enrichment = (candidate_with_term / length(candidate_locus_ids)) / 
                           (total_background_with_term / length(go_annotations$locus_id)),
          p_value = test_result$p.value,
          stringsAsFactors = FALSE
        ))
      }
    }
    
    if (nrow(results) > 0) {
      # Apply FDR correction
      results$p_adjusted <- p.adjust(results$p_value, method = "fdr")
      results <- results[order(results$p_value), ]
      
      cat("Tests performed:", nrow(results), "\n")
      cat("Significant at p < 0.05:", sum(results$p_value < 0.05), "\n")
      cat("Significant at FDR < 0.1:", sum(results$p_adjusted < 0.1), "\n")
      
      if (sum(results$p_adjusted < 0.1) > 0) {
        cat("Significant results (FDR < 0.1):\n")
        sig_results <- results[results$p_adjusted < 0.1, ]
        for (i in 1:nrow(sig_results)) {
          cat(sprintf("  %s: p=%.4f, FDR=%.4f, fold=%.2f\n",
                     sig_results$go_id[i], sig_results$p_value[i], 
                     sig_results$p_adjusted[i], sig_results$fold_enrichment[i]))
        }
      }
    }
    
    return(results)
  }
  
  # Test different subsets
  # 1. First 20 terms (like manual test)
  first_20 <- bp_terms_sized[1:min(20, length(bp_terms_sized))]
  results_20 <- test_subset(first_20, "First 20 terms (Manual Test Approach)")
  
  # 2. All terms (like new function)  
  results_all <- test_subset(bp_terms_sized, "All terms (New Function Approach)")
  
  # 3. Show the effect of multiple testing correction
  cat("\n=== Multiple Testing Effect Analysis ===\n")
  if (nrow(results_20) > 0 && nrow(results_all) > 0) {
    # Find common terms tested in both approaches
    common_terms <- intersect(results_20$go_id, results_all$go_id)
    
    if (length(common_terms) > 0) {
      cat("Common terms tested in both approaches:", length(common_terms), "\n")
      
      for (term in common_terms[1:min(5, length(common_terms))]) {
        p_20 <- results_20[results_20$go_id == term, "p_value"]
        fdr_20 <- results_20[results_20$go_id == term, "p_adjusted"]
        p_all <- results_all[results_all$go_id == term, "p_value"]
        fdr_all <- results_all[results_all$go_id == term, "p_adjusted"]
        
        cat(sprintf("  %s: p=%.4f (same), FDR_20terms=%.4f vs FDR_all=%.4f\n", 
                   term, p_20, fdr_20, fdr_all))
      }
    }
  }
  
  cat("\nConclusion: Testing more terms makes FDR correction more stringent,\n")
  cat("reducing the number of terms that reach significance.\n")
  
  return(list(
    results_20_terms = results_20,
    results_all_terms = results_all
  ))
}

# Usage:
# diagnosis <- diagnose_multiple_testing_effect(annotations, candidates$bed_file, "BP")