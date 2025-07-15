#!/usr/bin/env Rscript

# Check the actual database schema
library(funseqR)

cat("=== Checking Database Schema ===\n")

# Connect to database
con <- connect_funseq_db('funseq_project.db')

# Check if ora_analyses table exists and get its schema
tables <- DBI::dbListTables(con)
cat("Available tables:", paste(tables, collapse = ", "), "\n")

if ("ora_analyses" %in% tables) {
  cat("\nora_analyses table schema:\n")
  schema_info <- DBI::dbGetQuery(con, "PRAGMA table_info(ora_analyses)")
  print(schema_info)
  
  cat("\nSample query test (without parameters):\n")
  tryCatch({
    sample_data <- DBI::dbGetQuery(con, "SELECT * FROM ora_analyses LIMIT 1")
    cat("Query successful, rows returned:", nrow(sample_data), "\n")
  }, error = function(e) {
    cat("Query failed:", e$message, "\n")
  })
  
} else {
  cat("ora_analyses table does not exist!\n")
}

# Check blast_parameters table
if ("blast_parameters" %in% tables) {
  cat("\nblast_parameters table:\n")
  blast_data <- DBI::dbGetQuery(con, "SELECT blast_param_id FROM blast_parameters LIMIT 5")
  cat("blast_param_ids available:", paste(blast_data$blast_param_id, collapse = ", "), "\n")
} else {
  cat("blast_parameters table does not exist!\n")
}

# Close connection
close_funseq_db(con)