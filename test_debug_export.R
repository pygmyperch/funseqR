#!/usr/bin/env Rscript

# Test the flanking sequences export with debug output

library(funseqR)

cat("=== TESTING FLANKING EXPORT WITH DEBUG ===\n")

# Connect to source database
con <- connect_funseq_db('simple_funseq_project.db')

# Clean up any previous test database
test_db <- "test_debug_flanking.db"
if (file.exists(test_db)) {
  file.remove(test_db)
  cat("Removed existing test database\n")
}

tryCatch({
  # Test export with debug output enabled
  new_con <- export_project_data(
    con,
    test_db,
    export_candidates = FALSE,
    export_flanking = TRUE,
    export_blast = FALSE,
    export_annotations = FALSE,
    verbose = TRUE
  )
  
  if (!is.null(new_con)) {
    cat("✅ Export successful!\n")
    close_funseq_db(new_con)
  } else {
    cat("❌ Export failed\n")
  }
  
}, error = function(e) {
  cat("❌ Export error:\n")
  cat("Error message:", e$message, "\n")
  cat("Full error:\n")
  print(e)
})

# Close source connection
close_funseq_db(con)

cat("\n=== DEBUG TEST COMPLETE ===\n")