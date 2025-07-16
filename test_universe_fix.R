#!/usr/bin/env Rscript

# Test the ORA universe fix
# Compare results before/after fixing the clusterProfiler universe parameter

library(funseqR)

cat("=== TESTING ORA UNIVERSE FIX ===\n")
cat("Testing that fixed clusterProfiler universe gives consistent results\n\n")

# Connect to database
con <- connect_funseq_db('simple_funseq_project.db')

cat("1. Testing fixed clusterProfiler method...\n")
tryCatch({
  # Run ORA with fixed clusterProfiler method
  ora_results_fixed <- ora(con, 
                          annotation_type = "GO",
                          significance_threshold = 0.1,
                          method = "clusterprofiler",  # Uses fixed universe
                          store_results = FALSE,
                          create_plots = FALSE,
                          verbose = TRUE)
  
  cat("✓ Fixed clusterProfiler method completed successfully\n")
  
  # Extract significant terms
  if ("GO" %in% names(ora_results_fixed$enrichment_results)) {
    go_results <- ora_results_fixed$enrichment_results$GO
    for (ontology in c("BP", "MF", "CC")) {
      if (ontology %in% names(go_results)) {
        results <- go_results[[ontology]]
        if (nrow(results) > 0) {
          sig_terms <- results[results$p_adjusted < 0.1, ]
          cat("  ", ontology, "significant terms (FDR < 0.1):", nrow(sig_terms), "\n")
          if (nrow(sig_terms) > 0) {
            for (i in 1:min(3, nrow(sig_terms))) {
              term <- sig_terms[i, ]
              cat("    -", term$go_id, ":", term$go_term, "(p_adj =", 
                  format(term$p_adjusted, digits = 4), ")\n")
            }
          }
        }
      }
    }
  }
  
}, error = function(e) {
  cat("❌ Error with fixed clusterProfiler method:", e$message, "\n")
})

cat("\n2. Testing legacy method for comparison...\n")
tryCatch({
  # Run ORA with legacy method for comparison
  ora_results_legacy <- ora(con, 
                           annotation_type = "GO",
                           significance_threshold = 0.1,
                           method = "legacy",  # Uses hypergeometric test
                           store_results = FALSE,
                           create_plots = FALSE,
                           verbose = TRUE)
  
  cat("✓ Legacy method completed successfully\n")
  
  # Extract significant terms
  if ("GO" %in% names(ora_results_legacy$enrichment_results)) {
    go_results <- ora_results_legacy$enrichment_results$GO
    for (ontology in c("BP", "MF", "CC")) {
      if (ontology %in% names(go_results)) {
        results <- go_results[[ontology]]
        if (nrow(results) > 0) {
          sig_terms <- results[results$p_adjusted < 0.1, ]
          cat("  ", ontology, "significant terms (FDR < 0.1):", nrow(sig_terms), "\n")
          if (nrow(sig_terms) > 0) {
            for (i in 1:min(3, nrow(sig_terms))) {
              term <- sig_terms[i, ]
              cat("    -", term$go_id, ":", term$go_term, "(p_adj =", 
                  format(term$p_adjusted, digits = 4), ")\n")
            }
          }
        }
      }
    }
  }
  
}, error = function(e) {
  cat("❌ Error with legacy method:", e$message, "\n")
})

cat("\n3. Comparing methods...\n")

# Compare specific GO terms that should appear in both
if (exists("ora_results_fixed") && exists("ora_results_legacy")) {
  
  # Look for GO:0036342 (post-anal tail morphogenesis) in both
  cat("Checking GO:0036342 (post-anal tail morphogenesis):\n")
  
  # Check in fixed results
  fixed_go <- ora_results_fixed$enrichment_results$GO$BP
  if (!is.null(fixed_go) && nrow(fixed_go) > 0) {
    go_term <- fixed_go[fixed_go$go_id == "GO:0036342", ]
    if (nrow(go_term) > 0) {
      cat("  Fixed clusterProfiler: p_adj =", format(go_term$p_adjusted[1], digits = 6), 
          ", fold_enrichment =", format(go_term$fold_enrichment[1], digits = 3), "\n")
    } else {
      cat("  Fixed clusterProfiler: GO:0036342 not found\n")
    }
  }
  
  # Check in legacy results
  legacy_go <- ora_results_legacy$enrichment_results$GO$BP
  if (!is.null(legacy_go) && nrow(legacy_go) > 0) {
    go_term <- legacy_go[legacy_go$go_id == "GO:0036342", ]
    if (nrow(go_term) > 0) {
      cat("  Legacy method: p_adj =", format(go_term$p_adjusted[1], digits = 6),
          ", fold_enrichment =", format(go_term$fold_enrichment[1], digits = 3), "\n")
    } else {
      cat("  Legacy method: GO:0036342 not found\n")
    }
  }
  
  # Summary comparison
  cat("\nOverall comparison:\n")
  
  # Count significant terms in each method
  fixed_sig_count <- 0
  legacy_sig_count <- 0
  
  if (!is.null(ora_results_fixed$enrichment_results$GO)) {
    for (ont in names(ora_results_fixed$enrichment_results$GO)) {
      results <- ora_results_fixed$enrichment_results$GO[[ont]]
      if (!is.null(results) && nrow(results) > 0) {
        fixed_sig_count <- fixed_sig_count + sum(results$p_adjusted < 0.1, na.rm = TRUE)
      }
    }
  }
  
  if (!is.null(ora_results_legacy$enrichment_results$GO)) {
    for (ont in names(ora_results_legacy$enrichment_results$GO)) {
      results <- ora_results_legacy$enrichment_results$GO[[ont]]
      if (!is.null(results) && nrow(results) > 0) {
        legacy_sig_count <- legacy_sig_count + sum(results$p_adjusted < 0.1, na.rm = TRUE)
      }
    }
  }
  
  cat("  Fixed clusterProfiler significant terms:", fixed_sig_count, "\n")
  cat("  Legacy method significant terms:", legacy_sig_count, "\n")
  
  if (abs(fixed_sig_count - legacy_sig_count) <= 2) {
    cat("  ✅ Results are now consistent between methods!\n")
  } else {
    cat("  ⚠️  Results still differ - may need further investigation\n")
  }
}

cat("\n=== TEST COMPLETE ===\n")
cat("The universe fix should now provide statistically correct ORA results\n")
cat("that are consistent with the legacy hypergeometric method.\n")

# Close connection
close_funseq_db(con)