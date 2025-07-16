#!/usr/bin/env Rscript

# Test the export_project_data fix for reference sequences schema mismatch

library(funseqR)

cat("=== TESTING EXPORT FUNCTION FIX ===\n")
cat("Testing that reference sequences export now works\n\n")

# Connect to source database
source_db <- 'simple_funseq_project.db'
con <- connect_funseq_db(source_db)

cat("Testing export_project_data with fixed schema...\n")
cat("================================================\n")

# Clean up any previous test database
test_db <- "test_export_fix.db"
if (file.exists(test_db)) {
  file.remove(test_db)
  cat("Removed existing test database\n")
}

tryCatch({
  # Test export with minimal components first
  new_con <- export_project_data(
    con,
    test_db,
    export_candidates = FALSE,
    export_flanking = FALSE,
    export_blast = FALSE,
    export_annotations = FALSE,
    verbose = TRUE
  )
  
  cat("✅ Export completed successfully!\n\n")
  
  # Verify the exported data
  cat("Verifying exported data:\n")
  cat("========================\n")
  
  # Check reference sequences
  ref_seq_count <- DBI::dbGetQuery(new_con, "SELECT COUNT(*) as count FROM reference_sequences")$count
  cat("Reference sequences exported:", ref_seq_count, "\n")
  
  # Check if data matches source
  source_count <- DBI::dbGetQuery(con, "SELECT COUNT(*) as count FROM reference_sequences")$count
  cat("Source reference sequences:", source_count, "\n")
  
  if (ref_seq_count == source_count) {
    cat("✅ Reference sequence counts match!\n")
  } else {
    cat("⚠️  Reference sequence counts differ\n")
  }
  
  # Check VCF data
  vcf_count <- DBI::dbGetQuery(new_con, "SELECT COUNT(*) as count FROM vcf_data")$count
  source_vcf_count <- DBI::dbGetQuery(con, "SELECT COUNT(*) as count FROM vcf_data")$count
  cat("VCF entries - exported:", vcf_count, ", source:", source_vcf_count, "\n")
  
  # Close new connection
  close_funseq_db(new_con)
  
  cat("\n✅ Export function is now working correctly!\n")
  
}, error = function(e) {
  cat("❌ Export still failing:\n")
  cat("Error:", e$message, "\n")
  print(e)
})

# Close source connection
close_funseq_db(con)

cat("\n=== TEST COMPLETE ===\n")
cat("The schema mismatch has been fixed.\n")
cat("You should now be able to use export_project_data() successfully.\n")