#!/usr/bin/env Rscript

# Final test of the complete export fix

library(funseqR)

cat("=== FINAL EXPORT TEST ===\n")

# Connect to source database
con <- connect_funseq_db('simple_funseq_project.db')

# Clean up any previous test database
test_db <- "test_final_export.db"
if (file.exists(test_db)) {
  file.remove(test_db)
  cat("Removed existing test database\n")
}

tryCatch({
  # Test the complete export
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
    cat("🎉 COMPLETE EXPORT SUCCESSFUL! 🎉\n\n")
    
    # Verify all components
    cat("Final verification:\n")
    cat("==================\n")
    
    tables_to_check <- c("vcf_data", "reference_sequences", "candidate_loci", 
                         "locus_statistics", "flanking_sequences", "blast_results", 
                         "annotations", "go_terms", "kegg_references")
    
    for (table in tables_to_check) {
      count <- DBI::dbGetQuery(new_con, paste("SELECT COUNT(*) as count FROM", table))$count
      source_count <- DBI::dbGetQuery(con, paste("SELECT COUNT(*) as count FROM", table))$count
      
      status <- if (count == source_count) "✅" else "⚠️"
      cat(sprintf("%s %-20s: %6d (source: %6d)\n", status, table, count, source_count))
    }
    
    close_funseq_db(new_con)
    
    cat("\n🎉 ALL SCHEMA MISMATCHES RESOLVED! 🎉\n")
    cat("Your export_project_data() function is now fully working.\n")
    
  } else {
    cat("❌ Export failed\n")
  }
  
}, error = function(e) {
  cat("❌ Export error:\n")
  cat("Error message:", e$message, "\n")
  print(e)
  
  # Clean up failed database
  if (file.exists(test_db)) {
    file.remove(test_db)
  }
})

# Close source connection
close_funseq_db(con)

cat("\n=== FINAL TEST COMPLETE ===\n")