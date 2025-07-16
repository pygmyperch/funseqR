#!/usr/bin/env Rscript

# Debug the exact export error location

library(funseqR)
library(DBI)

cat("=== DEBUGGING EXPORT ERROR LOCATION ===\n")

# Connect to source database
source_db <- 'simple_funseq_project.db'
con <- connect_funseq_db(source_db)

# Clean up any previous test database
test_db <- "debug_export.db"
if (file.exists(test_db)) {
  file.remove(test_db)
  cat("Removed existing test database\n")
}

cat("Testing minimal export to isolate the error...\n")
cat("==============================================\n")

tryCatch({
  # Test just VCF + reference data first (this worked before)
  cat("Step 1: Testing VCF + Reference export...\n")
  new_con <- export_project_data(
    con,
    test_db,
    export_candidates = FALSE,
    export_flanking = FALSE,
    export_blast = FALSE,
    export_annotations = FALSE,
    verbose = TRUE
  )
  
  if (!is.null(new_con)) {
    cat("✅ VCF + Reference export successful\n")
    close_funseq_db(new_con)
    file.remove(test_db)
    
    # Now test adding candidates
    cat("\nStep 2: Testing VCF + Reference + Candidates export...\n")
    new_con2 <- export_project_data(
      con,
      test_db,
      export_candidates = TRUE,
      export_flanking = FALSE,
      export_blast = FALSE,
      export_annotations = FALSE,
      verbose = TRUE
    )
    
    if (!is.null(new_con2)) {
      cat("✅ VCF + Reference + Candidates export successful\n")
      close_funseq_db(new_con2)
      file.remove(test_db)
      
      # Now test the problematic flanking sequences
      cat("\nStep 3: Testing VCF + Reference + Flanking export (the problem)...\n")
      new_con3 <- export_project_data(
        con,
        test_db,
        export_candidates = FALSE,
        export_flanking = TRUE,
        export_blast = FALSE,
        export_annotations = FALSE,
        verbose = TRUE
      )
      
      if (!is.null(new_con3)) {
        cat("✅ Flanking sequences export successful!\n")
        close_funseq_db(new_con3)
      } else {
        cat("❌ Flanking sequences export failed\n")
      }
    } else {
      cat("❌ Candidates export failed\n")
    }
  } else {
    cat("❌ Basic export failed\n")
  }
  
}, error = function(e) {
  cat("❌ Export error:\n")
  cat("Error message:", e$message, "\n")
  cat("Full error details:\n")
  print(e)
  
  # Clean up
  if (file.exists(test_db)) {
    file.remove(test_db)
  }
})

# Close source connection
close_funseq_db(con)

cat("\n=== ERROR ISOLATION TEST COMPLETE ===\n")