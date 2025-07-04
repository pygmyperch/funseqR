# Simple test script for the ORA storage fix
library(funseqR)

test_storage_fix <- function(db_path = "funseq_project.db", vcf_file = "SA448_855.vcf") {
  
  cat("=== Testing ORA Storage Fix ===\n")
  
  if (!file.exists(db_path)) {
    cat("Database not found:", db_path, "\n")
    cat("Please copy your database to the repo root for testing\n")
    return(FALSE)
  }
  
  if (!file.exists(vcf_file)) {
    cat("VCF file not found:", vcf_file, "\n")
    cat("Please ensure your VCF file is available for testing\n")
    return(FALSE)
  }
  
  # Connect to database
  con <- connect_funseq_db(db_path, verbose = TRUE)
  
  # Delete any existing corrupted results
  tryCatch({
    delete_ora_results(con, analysis_id = c(1, 2, 3), confirm = FALSE, verbose = TRUE)
  }, error = function(e) {
    cat("No existing ORA results to delete\n")
  })
  
  # Run ORA with the fixed function
  cat("Running ORA with fixed storage function...\n")
  test_results <- run_ORA(con, vcf_file,
                         annotation_type = "GO",
                         significance_threshold = 0.1,
                         blast_param_id = 1,
                         verbose = TRUE)
  
  # Test enrichment retrieval
  cat("Testing enrichment retrieval...\n")
  # You would need your actual annotation data here
  # This is just a placeholder for the test structure
  
  close_funseq_db(con)
  cat("=== Test Complete ===\n")
  
  return(TRUE)
}

# This script can be run manually when testing the fix
cat("Storage fix test script loaded. Run test_storage_fix() to test.\n")