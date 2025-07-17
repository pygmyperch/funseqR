#!/usr/bin/env Rscript

# Debug what the source database actually contains

library(funseqR)
library(DBI)

cat("=== DEBUGGING SOURCE DATABASE SCHEMA ===\n")

# Connect to source database
con <- connect_funseq_db('simple_funseq_project.db')

cat("Checking actual schema of flanking_sequences table in source database...\n")
cat("====================================================================\n")

# Get actual table structure
schema_info <- DBI::dbGetQuery(con, "PRAGMA table_info(flanking_sequences)")
cat("Actual columns in source flanking_sequences table:\n")
print(schema_info)

cat("\nFirst few records to see data structure:\n")
sample_data <- DBI::dbGetQuery(con, "SELECT * FROM flanking_sequences LIMIT 3")
if (nrow(sample_data) > 0) {
  cat("Column names:", paste(colnames(sample_data), collapse = ", "), "\n")
  cat("Sample data:\n")
  print(sample_data)
} else {
  cat("No data in flanking_sequences table\n")
}

cat("\nAll tables in source database:\n")
cat("==============================\n")
tables <- DBI::dbListTables(con)
print(tables)

cat("\nChecking key tables schemas:\n")
cat("============================\n")

for (table in c("vcf_data", "reference_sequences", "blast_results", "annotations")) {
  if (table %in% tables) {
    cat("\n", table, ":\n")
    schema <- DBI::dbGetQuery(con, paste("PRAGMA table_info(", table, ")"))
    cat("Columns:", paste(schema$name, collapse = ", "), "\n")
  }
}

# Close connection
close_funseq_db(con)

cat("\n=== DEBUG COMPLETE ===\n")
cat("This will show us what the source database actually looks like\n")
cat("so we can understand the schema mismatch.\n")