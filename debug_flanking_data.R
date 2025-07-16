#!/usr/bin/env Rscript

# Debug the actual flanking data being exported

library(funseqR)
library(DBI)

cat("=== DEBUGGING FLANKING EXPORT DATA STRUCTURE ===\n")

# Connect to source database
con <- connect_funseq_db('simple_funseq_project.db')

cat("Testing the flanking export query...\n")
cat("===================================\n")

# This is the exact query used in .export_flanking_data
flanking_data <- DBI::dbGetQuery(con, "
  SELECT fs.*, vd.file_id as old_file_id
  FROM flanking_sequences fs
  JOIN vcf_data vd ON fs.vcf_id = vd.vcf_id
  ORDER BY fs.flanking_id
  LIMIT 3
")

cat("Query result structure:\n")
print(str(flanking_data))

cat("\nColumn names:\n")
print(colnames(flanking_data))

cat("\nFirst record:\n")
if (nrow(flanking_data) > 0) {
  print(flanking_data[1, ])
  
  cat("\nChecking for required columns:\n")
  required_cols <- c("vcf_id", "sequence_id", "flank_size", "start_position", "end_position", "sequence", "seq_type", "seq_length")
  
  for (col in required_cols) {
    if (col %in% colnames(flanking_data)) {
      cat("✅", col, ":", flanking_data[1, col], "\n")
    } else {
      cat("❌", col, ": MISSING\n")
    }
  }
} else {
  cat("No flanking data found!\n")
}

# Close connection
close_funseq_db(con)

cat("\n=== DEBUG COMPLETE ===\n")