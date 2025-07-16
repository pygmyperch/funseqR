#!/usr/bin/env Rscript

# Test the flanking sequences export fix

library(funseqR)

cat("=== TESTING FLANKING SEQUENCES EXPORT FIX ===\n")

# Connect to source database
con <- connect_funseq_db('simple_funseq_project.db')

# Clean up any previous test database
test_db <- "test_flanking_fix.db"
if (file.exists(test_db)) {
  file.remove(test_db)
  cat("Removed existing test database\n")
}

tryCatch({
  # Test export with just flanking sequences
  new_con <- export_project_data(
    con,
    test_db,
    export_candidates = FALSE,
    export_flanking = TRUE,
    export_blast = FALSE,
    export_annotations = FALSE,
    verbose = TRUE
  )
  
  cat("✅ Flanking sequences export completed successfully!\n\n")
  
  # Verify the exported data
  cat("Verifying exported flanking sequences:\n")
  cat("=====================================\n")
  
  # Check flanking sequences count
  flanking_count <- DBI::dbGetQuery(new_con, "SELECT COUNT(*) as count FROM flanking_sequences")$count
  source_flanking_count <- DBI::dbGetQuery(con, "SELECT COUNT(*) as count FROM flanking_sequences")$count
  cat("Flanking sequences - exported:", flanking_count, ", source:", source_flanking_count, "\n")
  
  if (flanking_count == source_flanking_count) {
    cat("✅ All flanking sequences exported successfully!\n")
  } else {
    cat("⚠️  Some flanking sequences missing in export\n")
  }
  
  # Close new connection
  close_funseq_db(new_con)
  
}, error = function(e) {
  cat("❌ Flanking sequences export still failing:\n")
  cat("Error:", e$message, "\n")
  print(e)
})

# Close source connection
close_funseq_db(con)

cat("\n=== TEST COMPLETE ===\n")