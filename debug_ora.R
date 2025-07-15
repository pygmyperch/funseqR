#!/usr/bin/env Rscript

# Debug the ORA storage issue
library(funseqR)

cat("=== Debugging ORA Storage Issue ===\n")

# Connect to database
con <- connect_funseq_db('funseq_project.db')

# Run ORA with minimal settings to trigger the error
cat("\nRunning ORA to trigger storage error...\n")
tryCatch({
  ora_result <- ora(con, 
                   candidate_vcf_file = "stored",
                   annotation_type = "GO",
                   ontologies = c("BP"),  # Just test BP to isolate the issue
                   significance_threshold = 0.1,
                   store_results = TRUE,   # This should trigger the error
                   create_plots = FALSE,
                   verbose = TRUE)
  cat("SUCCESS: ORA completed\n")
}, error = function(e) {
  cat("ERROR in ORA:", e$message, "\n")
  cat("ERROR details:\n")
  print(e)
})

# Close connection
close_funseq_db(con)