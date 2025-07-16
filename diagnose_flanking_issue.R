#!/usr/bin/env Rscript

# Diagnose the flanking sequences export issue

library(DBI)
library(RSQLite)

cat("=== DIAGNOSING FLANKING SEQUENCES EXPORT ISSUE ===\n")

# Connect to the database
con <- dbConnect(SQLite(), 'simple_funseq_project.db')

# Check the table schema
cat("1. Checking flanking_sequences table schema:\n")
table_info <- dbGetQuery(con, "PRAGMA table_info(flanking_sequences)")
print(table_info)

cat("\n2. Checking actual data structure:\n")
count_result <- dbGetQuery(con, "SELECT COUNT(*) as count FROM flanking_sequences")
cat("Number of rows in flanking_sequences:", count_result$count, "\n")

if (count_result$count > 0) {
  # Get first row to see actual data
  sample_data <- dbGetQuery(con, "SELECT * FROM flanking_sequences LIMIT 1")
  cat("\nSample data:\n")
  print(sample_data)
  
  # Check data types
  cat("\nData types:\n")
  for (col in colnames(sample_data)) {
    cat(col, ":", class(sample_data[[col]]), "\n")
  }
}

cat("\n3. Analyzing the INSERT statement mismatch:\n")
cat("INSERT statement expects 8 parameters for columns:\n")
expected_cols <- c("vcf_id", "sequence_id", "flank_size", "start_position", "end_position", "sequence", "seq_type", "seq_length")
cat(paste(expected_cols, collapse = ", "), "\n")

actual_cols <- table_info$name
cat("\nActual table columns:\n")
cat(paste(actual_cols, collapse = ", "), "\n")

# Check if the INSERT statement is missing the PRIMARY KEY
pk_col <- table_info$name[table_info$pk == 1]
cat("\nPrimary key column:", pk_col, "\n")

if (!(pk_col %in% expected_cols)) {
  cat("❌ PRIMARY KEY column", pk_col, "is missing from INSERT statement!\n")
  cat("This is likely the source of the parameter mismatch.\n")
}

# Check for any NULL constraints that might affect the INSERT
cat("\n4. Column constraints:\n")
for (i in 1:nrow(table_info)) {
  col_info <- table_info[i, ]
  if (col_info$notnull == 1) {
    cat("Column", col_info$name, "has NOT NULL constraint\n")
  }
}

# Test the specific query that's causing the issue
cat("\n5. Testing the export query:\n")
tryCatch({
  flanking_data <- dbGetQuery(con, "
    SELECT fs.*, vd.file_id as old_file_id
    FROM flanking_sequences fs
    JOIN vcf_data vd ON fs.vcf_id = vd.vcf_id
    ORDER BY fs.flanking_id
    LIMIT 1
  ")
  
  if (nrow(flanking_data) > 0) {
    cat("Export query columns:\n")
    print(colnames(flanking_data))
    
    # Check if we're trying to insert a record that includes the primary key
    if ("flanking_id" %in% colnames(flanking_data)) {
      cat("❌ flanking_id is present in export data - this will cause issues\n")
      cat("The INSERT statement doesn't include flanking_id in the column list\n")
    }
    
    # Simulate the record construction
    fs_record <- flanking_data[1, ]
    
    # Check parameter lengths
    cat("\nParameter analysis:\n")
    params <- list(
      fs_record$vcf_id,
      fs_record$sequence_id,
      fs_record$flank_size,
      fs_record$start_position,
      fs_record$end_position,
      fs_record$sequence,
      fs_record$seq_type,
      fs_record$seq_length
    )
    
    for (i in 1:length(params)) {
      param <- params[[i]]
      cat("Parameter", i, ":", param, "- length:", length(param), "- class:", class(param), "\n")
      
      if (length(param) != 1) {
        cat("❌ Parameter", i, "has length", length(param), "instead of 1!\n")
      }
    }
  }
  
}, error = function(e) {
  cat("❌ Error in export query:", e$message, "\n")
})

dbDisconnect(con)
cat("\n=== DIAGNOSIS COMPLETE ===\n")