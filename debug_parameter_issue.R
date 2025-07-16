#!/usr/bin/env Rscript

# Debug the parameter length issue in flanking sequences

library(funseqR)
library(DBI)

cat("=== DEBUGGING PARAMETER LENGTH ISSUE ===\n")

# Connect to source database
con <- connect_funseq_db('simple_funseq_project.db')

cat("Analyzing the flanking_data structure being passed to batch insert...\n")
cat("===================================================================\n")

# Simulate the exact data preparation that happens in .export_flanking_data
flanking_data <- DBI::dbGetQuery(con, "
  SELECT fs.*, vd.file_id as old_file_id
  FROM flanking_sequences fs
  JOIN vcf_data vd ON fs.vcf_id = vd.vcf_id
  ORDER BY fs.flanking_id
  LIMIT 5
")

cat("Source flanking_data columns:\n")
print(colnames(flanking_data))

if (nrow(flanking_data) > 0) {
  cat("\nFirst record structure:\n")
  print(str(flanking_data[1, ]))
  
  # Simulate the export_record creation
  fs_record <- flanking_data[1, ]
  
  cat("\nChecking individual field values:\n")
  cat("vcf_id:", fs_record$vcf_id, "- class:", class(fs_record$vcf_id), "- length:", length(fs_record$vcf_id), "\n")
  
  # Check if these columns exist in the source data
  expected_cols <- c("sequence_id", "flank_size", "start_position", "end_position", "sequence", "seq_type", "seq_length")
  
  for (col in expected_cols) {
    if (col %in% colnames(fs_record)) {
      value <- fs_record[[col]]
      cat(col, ":", value, "- class:", class(value), "- length:", length(value), "- is.na:", is.na(value), "\n")
    } else {
      cat("❌", col, ": COLUMN MISSING FROM SOURCE DATA\n")
    }
  }
  
  # Check what would be the export_record structure
  cat("\nTesting export_record construction...\n")
  
  # Check if we can access the fields without error
  tryCatch({
    if ("sequence_id" %in% colnames(fs_record)) {
      export_record <- data.frame(
        vcf_id = 1,  # dummy value
        sequence_id = fs_record$sequence_id,
        flank_size = if("flank_size" %in% colnames(fs_record)) fs_record$flank_size else NA,
        start_position = if("start_position" %in% colnames(fs_record)) fs_record$start_position else NA,
        end_position = if("end_position" %in% colnames(fs_record)) fs_record$end_position else NA,
        sequence = if("sequence" %in% colnames(fs_record)) fs_record$sequence else NA,
        seq_type = if("seq_type" %in% colnames(fs_record)) fs_record$seq_type else "raw",
        seq_length = if("seq_length" %in% colnames(fs_record)) fs_record$seq_length else NA
      )
      
      cat("✅ Export record created successfully:\n")
      print(str(export_record))
      
      cat("\nChecking for NULL or NA values that could cause parameter issues:\n")
      for (i in 1:ncol(export_record)) {
        col_name <- colnames(export_record)[i]
        value <- export_record[1, i]
        cat("Parameter", i, "(", col_name, "):", value, "- is.na:", is.na(value), "- is.null:", is.null(value), "\n")
      }
      
    } else {
      cat("❌ sequence_id column missing - cannot create export record\n")
    }
    
  }, error = function(e) {
    cat("❌ Error creating export record:", e$message, "\n")
  })
  
} else {
  cat("No flanking data found!\n")
}

# Close connection
close_funseq_db(con)

cat("\n=== PARAMETER DEBUG COMPLETE ===\n")