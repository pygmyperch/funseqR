#!/usr/bin/env Rscript

# Quick check of actual flanking_sequences columns

library(funseqR)
library(DBI)

cat("=== QUICK COLUMN CHECK ===\n")

# Connect to source database
con <- connect_funseq_db('simple_funseq_project.db')

# Check table structure
cat("Actual flanking_sequences table structure:\n")
schema_info <- DBI::dbGetQuery(con, "PRAGMA table_info(flanking_sequences)")
print(schema_info)

cat("\nSample data (first row):\n")
sample_data <- DBI::dbGetQuery(con, "SELECT * FROM flanking_sequences LIMIT 1")
if (nrow(sample_data) > 0) {
  print(sample_data)
  cat("\nColumn names:\n")
  print(colnames(sample_data))
} else {
  cat("No data found!\n")
}

# Close connection
close_funseq_db(con)

cat("\n=== CHECK COMPLETE ===\n")