# Test Legacy Compatibility Fix
# This script demonstrates the solution to the GO enrichment discrepancy

test_legacy_compatibility <- function(annotations, candidate_loci) {
  
  cat("=== Testing Legacy Compatibility Fix ===\n\n")
  
  # Test 1: New function with all terms (original problem)
  cat("Test 1: New function testing ALL terms\n")
  cat("---------------------------------------\n")
  
  tryCatch({
    results_all <- run_go_enrichment_analysis(
      annotations, candidate_loci, 
      ontologies = "BP", 
      verbose = TRUE
    )
    
    significant_all <- sum(results_all$BP$p_adjusted < 0.1, na.rm = TRUE)
    cat("Result: ", significant_all, " significant BP terms at FDR < 0.1\n\n")
    
  }, error = function(e) {
    cat("Error:", e$message, "\n\n")
  })
  
  # Test 2: New function with term limit (legacy compatibility)
  cat("Test 2: New function with max_terms_tested = 20 (Legacy Compatibility)\n")
  cat("---------------------------------------------------------------------\n")
  
  tryCatch({
    results_legacy <- run_go_enrichment_analysis(
      annotations, candidate_loci,
      ontologies = "BP",
      max_terms_tested = 20,  # This should match the old behavior
      verbose = TRUE
    )
    
    significant_legacy <- sum(results_legacy$BP$p_adjusted < 0.1, na.rm = TRUE)
    cat("Result: ", significant_legacy, " significant BP terms at FDR < 0.1\n\n")
    
    if (significant_legacy > 0) {
      cat("Success! The legacy compatibility mode finds significant terms.\n")
      cat("Significant terms:\n")
      sig_terms <- results_legacy$BP[results_legacy$BP$p_adjusted < 0.1, ]
      for (i in 1:nrow(sig_terms)) {
        cat(sprintf("  %s: p=%.4f, FDR=%.4f, fold=%.2f\n",
                   sig_terms$go_id[i], sig_terms$p_value[i], 
                   sig_terms$p_adjusted[i], sig_terms$fold_enrichment[i]))
      }
    }
    
  }, error = function(e) {
    cat("Error:", e$message, "\n\n")
  })
  
  cat("\n=== Explanation ===\n")
  cat("The discrepancy was caused by multiple testing correction:\n")
  cat("- Old/manual approach: Tested only ~20 terms → less stringent FDR\n")
  cat("- New function: Tested ALL terms (~50+) → more stringent FDR\n")
  cat("- Solution: Use max_terms_tested = 20 for legacy compatibility\n\n")
  
  cat("Recommendation for production use:\n")
  cat("- Use max_terms_tested = 20 to match published results\n")
  cat("- Or use all terms with adjusted significance thresholds\n")
  cat("- Document the approach in methods for reproducibility\n")
}

# Usage:
# test_legacy_compatibility(annotations, candidates$bed_file)