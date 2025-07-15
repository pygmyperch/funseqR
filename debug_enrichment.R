#!/usr/bin/env Rscript

# Debug the enrichment results structure
library(funseqR)

cat("=== Debugging Enrichment Results Structure ===\n")

# Connect to database
con <- connect_funseq_db('funseq_project.db')

# Run just the enrichment analysis part without storage
cat("\nExtracting GO data...\n")
go_data <- extract_go_terms_for_enrichment(con, "stored", verbose = TRUE)

cat("\nPerforming GO enrichment (BP only)...\n")
bp_results <- perform_go_enrichment(go_data, "BP", min_genes = 5, max_genes = 500, 
                                  significance_threshold = 0.1, method = "clusterprofiler", 
                                  verbose = TRUE)

cat("\nBP Results structure:\n")
cat("Class:", class(bp_results), "\n")
cat("Dimensions:", nrow(bp_results), "x", ncol(bp_results), "\n")
cat("Column names:", paste(colnames(bp_results), collapse = ", "), "\n")

if (nrow(bp_results) > 0) {
  cat("\nFirst row values:\n")
  first_row <- bp_results[1, ]
  for (col in colnames(bp_results)) {
    cat("  ", col, ":", first_row[[col]], "\n")
  }
}

# Close connection
close_funseq_db(con)