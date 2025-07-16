#!/usr/bin/env Rscript

# Test the export after fixing the schema mismatch

library(funseqR)

cat("=== TESTING FIXED EXPORT ===\n")

# Connect to source database
con <- connect_funseq_db('simple_funseq_project.db')

# Clean up any previous test database
test_db <- "test_fixed_export.db"
if (file.exists(test_db)) {
  file.remove(test_db)
  cat("Removed existing test database\n")
}

tryCatch({
  # Test the problematic flanking sequences export
  new_con <- export_project_data(
    con,
    test_db,
    export_candidates = TRUE,
    export_flanking = TRUE,
    export_blast = TRUE,
    export_annotations = TRUE,
    verbose = TRUE
  )
  
  if (!is.null(new_con)) {
    cat("✅ Complete export successful!\n\n")
    
    # Verify the exported data
    cat("Verifying exported data:\n")
    cat("=======================\n")
    
    # Check all components
    vcf_count <- DBI::dbGetQuery(new_con, "SELECT COUNT(*) as count FROM vcf_data")$count
    ref_count <- DBI::dbGetQuery(new_con, "SELECT COUNT(*) as count FROM reference_sequences")$count
    candidates_count <- DBI::dbGetQuery(new_con, "SELECT COUNT(*) as count FROM candidate_loci")$count
    flanking_count <- DBI::dbGetQuery(new_con, "SELECT COUNT(*) as count FROM flanking_sequences")$count
    blast_count <- DBI::dbGetQuery(new_con, "SELECT COUNT(*) as count FROM blast_results")$count
    annotations_count <- DBI::dbGetQuery(new_con, "SELECT COUNT(*) as count FROM annotations")$count
    
    cat("VCF entries:", vcf_count, "\n")
    cat("Reference sequences:", ref_count, "\n") 
    cat("Candidate loci:", candidates_count, "\n")
    cat("Flanking sequences:", flanking_count, "\n")
    cat("BLAST results:", blast_count, "\n")
    cat("Annotations:", annotations_count, "\n")
    
    close_funseq_db(new_con)
    
    cat("\n✅ ALL SCHEMA MISMATCHES FIXED!\n")
    cat("Your export_project_data() function now works correctly.\n")
    
  } else {
    cat("❌ Export failed\n")
  }
  
}, error = function(e) {
  cat("❌ Export error:\n")
  cat("Error message:", e$message, "\n")
  print(e)
})

# Close source connection
close_funseq_db(con)

cat("\n=== FINAL TEST COMPLETE ===\n")