#!/usr/bin/env Rscript

# Debug the flanking sequences schema mismatch

library(funseqR)
library(DBI)

cat("=== DEBUGGING FLANKING SEQUENCES SCHEMA ===\n")

# Connect to source database
con <- connect_funseq_db('simple_funseq_project.db')

cat("Checking actual flanking_sequences table structure...\n")
cat("================================================\n")

# Get table schema
schema_info <- DBI::dbGetQuery(con, "PRAGMA table_info(flanking_sequences)")
print(schema_info)

cat("\nFirst few records from flanking_sequences:\n")
cat("==========================================\n")
sample_data <- DBI::dbGetQuery(con, "SELECT * FROM flanking_sequences LIMIT 3")
print(sample_data)

cat("\nColumn names in flanking_sequences:\n")
cat("==================================\n")
print(colnames(sample_data))

cat("\nChecking if columns match schema expectations:\n")
cat("==============================================\n")
expected_cols <- c("flanking_id", "vcf_id", "sequence_id", "flank_size", "start_position", "end_position", "sequence", "seq_type", "seq_length")
actual_cols <- colnames(sample_data)

cat("Expected columns:", paste(expected_cols, collapse = ", "), "\n")
cat("Actual columns:  ", paste(actual_cols, collapse = ", "), "\n")

missing_cols <- setdiff(expected_cols, actual_cols)
extra_cols <- setdiff(actual_cols, expected_cols)

if (length(missing_cols) > 0) {
  cat("❌ Missing columns:", paste(missing_cols, collapse = ", "), "\n")
}
if (length(extra_cols) > 0) {
  cat("⚠️  Extra columns:", paste(extra_cols, collapse = ", "), "\n")
}

if (length(missing_cols) == 0 && length(extra_cols) == 0) {
  cat("✅ Column structure matches expected schema\n")
}

# Close connection
close_funseq_db(con)

cat("\n=== DEBUG COMPLETE ===\n")