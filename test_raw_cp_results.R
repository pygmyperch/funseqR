#!/usr/bin/env Rscript

# Test access to full clusterProfiler results from updated ora() function
library(funseqR)

cat("=== Testing Raw clusterProfiler Results Access ===\n")

# Connect to database
con <- connect_funseq_db('funseq_project.db')

cat("\n1. Running ORA with updated function...\n")
tryCatch({
  # Run ORA analysis 
  ora_results <- ora(con, 
                    annotation_type = "GO",
                    significance_threshold = 0.1,
                    store_results = FALSE,  # Don't store for this test
                    create_plots = FALSE,   # Don't create plots for this test
                    verbose = TRUE)
  
  cat("SUCCESS: ORA function completed\n")
  cat("  - Status:", ora_results$status, "\n")
  
  # Check if raw_clusterprofiler_results is available
  if ("raw_clusterprofiler_results" %in% names(ora_results)) {
    cat("  - Raw clusterProfiler results available: YES\n")
    
    # Check GO results structure
    if ("GO" %in% names(ora_results$raw_clusterprofiler_results)) {
      go_raw <- ora_results$raw_clusterprofiler_results$GO
      cat("  - GO raw results available: YES\n")
      
      # Check each ontology
      for (ontology in c("BP", "MF", "CC")) {
        if (ontology %in% names(go_raw)) {
          raw_result <- go_raw[[ontology]]
          cat("    - ", ontology, " raw object class:", class(raw_result)[1], "\n")
          
          # Try to access clusterProfiler-specific data
          if (!is.null(raw_result) && "enrichResult" %in% class(raw_result)) {
            cat("    - ", ontology, " is enrichResult object: YES\n")
            if (nrow(raw_result@result) > 0) {
              cat("    - ", ontology, " enriched terms:", nrow(raw_result@result), "\n")
              # Show that we can access full clusterProfiler data
              result_df <- as.data.frame(raw_result)
              cat("    - ", ontology, " columns available:", ncol(result_df), "\n")
              cat("    - ", ontology, " column names:", paste(names(result_df)[1:min(5, ncol(result_df))], collapse = ", "), "...\n")
            } else {
              cat("    - ", ontology, " has no enriched terms\n")
            }
          } else {
            cat("    - ", ontology, " raw object:", if(is.null(raw_result)) "NULL" else "Not enrichResult", "\n")
          }
        }
      }
    }
    
    # Compare with funseqR results
    if ("enrichment_results" %in% names(ora_results) && "GO" %in% names(ora_results$enrichment_results)) {
      funseqr_go <- ora_results$enrichment_results$GO
      cat("\n2. Comparing with funseqR formatted results:\n")
      
      for (ontology in c("BP", "MF", "CC")) {
        if (ontology %in% names(funseqr_go) && ontology %in% names(go_raw)) {
          funseqr_terms <- nrow(funseqr_go[[ontology]])
          raw_terms <- if (!is.null(go_raw[[ontology]])) nrow(go_raw[[ontology]]@result) else 0
          cat("  - ", ontology, " terms - funseqR:", funseqr_terms, "raw:", raw_terms, "\n")
        }
      }
    }
    
  } else {
    cat("  - Raw clusterProfiler results available: NO\n")
  }
  
}, error = function(e) {
  cat("ERROR in ORA function:\n")
  cat("  Message:", e$message, "\n")
  cat("  Class:", class(e), "\n")
  print(e)
})

cat("\n=== Test Complete ===\n")

# Close connection
close_funseq_db(con)