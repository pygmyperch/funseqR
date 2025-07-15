#!/usr/bin/env Rscript

# Test the final ORA fix
library(funseqR)

cat("=== Testing Final ORA Fix ===\n")

# Connect to database
con <- connect_funseq_db('funseq_project.db')

# Run ORA with the fixed storage function
cat("\nRunning ORA with fixed storage...\n")
tryCatch({
  ora_result <- ora(con, 
                   candidate_vcf_file = "stored",
                   annotation_type = "GO",
                   ontologies = c("BP"),  # Just test BP first
                   significance_threshold = 0.1,
                   store_results = TRUE,   # Test storage with the fix
                   create_plots = FALSE,
                   verbose = TRUE)
  cat("SUCCESS: ORA completed with storage!\n")
  cat("  - Status:", ora_result$status, "\n")
  if (!is.null(ora_result$analysis_ids)) {
    cat("  - Analysis IDs:", paste(unlist(ora_result$analysis_ids), collapse = ", "), "\n")
  }
  if (!is.null(ora_result$enrichment_results$GO$BP)) {
    cat("  - BP enrichments found:", nrow(ora_result$enrichment_results$GO$BP), "\n")
  }
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  print(e)
})

# Close connection
close_funseq_db(con)