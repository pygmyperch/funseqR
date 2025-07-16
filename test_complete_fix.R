#!/usr/bin/env Rscript

# Test the complete ORA fix: universe + TERM2GENE mapping
# This should now give results consistent with the cp_enrich database

library(funseqR)

cat("=== TESTING COMPLETE ORA FIXES ===\n")
cat("Testing universe + TERM2GENE mapping fixes\n\n")

# Connect to database
con <- connect_funseq_db('simple_funseq_project.db')

cat("Running ORA with all fixes applied...\n")
cat("====================================\n")

tryCatch({
  # Run ORA with fixed clusterProfiler method
  ora_results <- ora(con, 
                    annotation_type = "GO",
                    significance_threshold = 0.1,
                    method = "clusterprofiler",  # Uses both fixes
                    store_results = FALSE,
                    create_plots = FALSE,
                    verbose = TRUE)
  
  cat("✅ Fixed clusterProfiler method completed successfully\n\n")
  
  # Analyze results in detail
  if ("GO" %in% names(ora_results$enrichment_results)) {
    go_results <- ora_results$enrichment_results$GO
    
    # Check BP results specifically
    if ("BP" %in% names(go_results)) {
      bp_results <- go_results[["BP"]]
      if (nrow(bp_results) > 0) {
        
        # Look for the terms that appeared in your original comparison
        target_terms <- c("GO:0036342", "GO:0030097")  # post-anal tail, hemopoiesis
        
        cat("Checking target GO terms:\n")
        cat("========================\n")
        
        for (term_id in target_terms) {
          term_row <- bp_results[bp_results$go_id == term_id, ]
          if (nrow(term_row) > 0) {
            cat("✓", term_id, "found:\n")
            cat("  - Term:", term_row$go_term[1], "\n")
            cat("  - p_value:", format(term_row$p_value[1], scientific = TRUE, digits = 4), "\n")
            cat("  - p_adjusted:", format(term_row$p_adjusted[1], scientific = TRUE, digits = 4), "\n")
            cat("  - Significant (FDR < 0.1):", term_row$p_adjusted[1] < 0.1, "\n")
            cat("  - Foreground count:", term_row$foreground_count[1], "\n")
            cat("  - Background count:", term_row$background_count[1], "\n")
            cat("  - Fold enrichment:", format(term_row$fold_enrichment[1], digits = 3), "\n")
            
            # Check gene IDs if available
            if ("gene_ids" %in% colnames(term_row) && !is.na(term_row$gene_ids[1])) {
              genes <- strsplit(term_row$gene_ids[1], "/")[[1]]
              cat("  - Genes:", paste(head(genes, 3), collapse = ", "), 
                  if(length(genes) > 3) paste("... (", length(genes), "total)") else "", "\n")
            }
            cat("\n")
          } else {
            cat("✗", term_id, "not found or not significant\n\n")
          }
        }
        
        # Show all significant BP terms
        sig_bp <- bp_results[bp_results$p_adjusted < 0.1, ]
        cat("All significant BP terms (FDR < 0.1):\n")
        cat("====================================\n")
        
        if (nrow(sig_bp) > 0) {
          for (i in 1:nrow(sig_bp)) {
            term <- sig_bp[i, ]
            cat(i, ". ", term$go_id, ": ", term$go_term, "\n", sep = "")
            cat("   p_adj = ", format(term$p_adjusted, scientific = TRUE, digits = 4), 
                ", fold = ", format(term$fold_enrichment, digits = 3), "\n")
          }
        } else {
          cat("No significant BP terms found\n")
        }
        
        cat("\nThis should now match the cp_enrich database results!\n")
        
      } else {
        cat("No BP results found\n")
      }
    }
  }
  
}, error = function(e) {
  cat("❌ Error with fixed clusterProfiler method:", e$message, "\n")
  cat("Full error:\n")
  print(e)
})

cat("\n=== EXPECTED RESULTS ===\n")
cat("If the fixes are correct, you should now see:\n")
cat("1. GO:0036342 (post-anal tail morphogenesis) as significant\n")
cat("2. GO:0030097 (hemopoiesis) as significant\n") 
cat("3. Similar p_adjusted values to the cp_enrich database\n")
cat("4. The same enriched loci as the cp_enrich version\n")

cat("\n=== TEST COMPLETE ===\n")

# Close connection
close_funseq_db(con)