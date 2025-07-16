#!/usr/bin/env Rscript

# Check the actual schema of the flanking_sequences table

library(DBI)
library(RSQLite)

# Connect to the database
con <- dbConnect(SQLite(), 'simple_funseq_project.db')

cat("=== CHECKING FLANKING_SEQUENCES TABLE SCHEMA ===\n")

# Get table info
table_info <- dbGetQuery(con, "PRAGMA table_info(flanking_sequences)")
cat("Table schema:\n")
print(table_info)

cat("\n=== CHECKING ACTUAL DATA STRUCTURE ===\n")

# Check if table exists and has data
count_query <- "SELECT COUNT(*) as count FROM flanking_sequences"
count_result <- dbGetQuery(con, count_query)
cat("Number of rows in flanking_sequences:", count_result$count, "\n")

if (count_result$count > 0) {
  # Get first few rows to see actual data
  sample_data <- dbGetQuery(con, "SELECT * FROM flanking_sequences LIMIT 3")
  cat("\nSample data columns:\n")
  print(colnames(sample_data))
  
  cat("\nSample data structure:\n")
  print(str(sample_data))
  
  cat("\nFirst row values:\n")
  if (nrow(sample_data) > 0) {
    print(sample_data[1, ])
  }
}

# Check what the batch insert expects vs what's in the table
cat("\n=== EXPECTED VS ACTUAL COLUMNS ===\n")
expected_cols <- c("vcf_id", "sequence_id", "flank_size", "start_position", "end_position", "sequence", "seq_type", "seq_length")
actual_cols <- table_info$name

cat("Expected columns for INSERT:", paste(expected_cols, collapse = ", "), "\n")
cat("Actual table columns:", paste(actual_cols, collapse = ", "), "\n")

missing_cols <- setdiff(expected_cols, actual_cols)
extra_cols <- setdiff(actual_cols, expected_cols)

if (length(missing_cols) > 0) {
  cat("❌ Missing columns:", paste(missing_cols, collapse = ", "), "\n")
}

if (length(extra_cols) > 0) {
  cat("➕ Extra columns:", paste(extra_cols, collapse = ", "), "\n")
}

# Check the INSERT statement parameter count
cat("\n=== INSERT STATEMENT ANALYSIS ===\n")
cat("INSERT statement expects 8 parameters (?, ?, ?, ?, ?, ?, ?, ?)\n")
cat("Corresponding to: vcf_id, sequence_id, flank_size, start_position, end_position, sequence, seq_type, seq_length\n")
cat("Table has", nrow(table_info), "columns\n")

# Check if flanking_id is auto-increment
primary_key <- table_info[table_info$pk == 1, "name"]
cat("Primary key column:", primary_key, "\n")

dbDisconnect(con)
cat("\n=== SCHEMA CHECK COMPLETE ===\n")