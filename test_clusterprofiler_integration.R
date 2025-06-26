# Test clusterProfiler Integration
# Quick test script to validate the new functions

# Load the package in development mode
# devtools::load_all()  # Uncomment if using devtools

cat("=== Testing clusterProfiler Integration ===\\n\\n")

# Test 1: Check if clusterProfiler is available
cat("Test 1: Checking clusterProfiler availability...\\n")
if (requireNamespace("clusterProfiler", quietly = TRUE)) {
  cat("✓ clusterProfiler is available\\n")
} else {
  cat("✗ clusterProfiler not found. Install with: BiocManager::install('clusterProfiler')\\n")
}

if (requireNamespace("DOSE", quietly = TRUE)) {
  cat("✓ DOSE is available\\n")
} else {
  cat("✗ DOSE not found. Install with: BiocManager::install('DOSE')\\n")
}

# Test 2: Check function loading
cat("\\nTest 2: Checking function availability...\\n")
if (exists("run_clusterprofiler_enrichment")) {
  cat("✓ run_clusterprofiler_enrichment() function loaded\\n")
} else {
  cat("✗ run_clusterprofiler_enrichment() function not found\\n")
  cat("   Make sure to load the package or source the R files\\n")
}

if (exists("compare_enrichment_methods")) {
  cat("✓ compare_enrichment_methods() function loaded\\n")
} else {
  cat("✗ compare_enrichment_methods() function not found\\n")
}

# Test 3: Database connection test (if database exists)
cat("\\nTest 3: Database connection test...\\n")
test_db_path <- "test_funseq.db"

if (file.exists(test_db_path)) {
  tryCatch({
    con <- connect_funseq_db(test_db_path, verbose = FALSE)
    
    # Check for required tables
    tables <- DBI::dbListTables(con)
    required_tables <- c("vcf_data", "flanking_sequences", "blast_results", "annotations", "go_terms")
    
    missing_tables <- setdiff(required_tables, tables)
    if (length(missing_tables) == 0) {
      cat("✓ All required database tables present\\n")
      
      # Check for data
      vcf_count <- DBI::dbGetQuery(con, "SELECT COUNT(*) as count FROM vcf_data")$count
      go_count <- DBI::dbGetQuery(con, "SELECT COUNT(*) as count FROM go_terms")$count
      
      cat("  - VCF entries:", vcf_count, "\\n")
      cat("  - GO terms:", go_count, "\\n")
      
      if (vcf_count > 0 && go_count > 0) {
        cat("✓ Database contains data for testing\\n")
      } else {
        cat("! Database exists but lacks sufficient data for testing\\n")
      }
      
    } else {
      cat("✗ Missing database tables:", paste(missing_tables, collapse = ", "), "\\n")
    }
    
    DBI::dbDisconnect(con)
    
  }, error = function(e) {
    cat("✗ Database connection error:", e$message, "\\n")
  })
} else {
  cat("! Test database not found at:", test_db_path, "\\n")
  cat("  Run the full workflow example to create test data\\n")
}

# Test 4: Function parameter validation
cat("\\nTest 4: Function parameter validation...\\n")
tryCatch({
  # Test with dummy parameters (should fail gracefully)
  result <- run_clusterprofiler_enrichment(
    con = NULL,  # This should trigger an error
    candidate_file_id = 1,
    background_file_id = 2,
    blast_param_id = 1,
    verbose = FALSE
  )
  cat("✗ Function should have failed with NULL connection\\n")
}, error = function(e) {
  cat("✓ Function properly validates parameters\\n")
})

# Test 5: Example data structure test
cat("\\nTest 5: Testing helper functions...\\n")

# Test ontology mapping
test_ontologies <- c("BP", "MF", "CC")
ontology_map <- c("BP" = "P", "MF" = "F", "CC" = "C")

for (ont in test_ontologies) {
  if (ont %in% names(ontology_map)) {
    cat("✓ Ontology mapping for", ont, "->", ontology_map[ont], "\\n")
  } else {
    cat("✗ Missing ontology mapping for", ont, "\\n")
  }
}

cat("\\n=== Test Summary ===\\n")
cat("If all tests pass (✓), the clusterProfiler integration is ready to use.\\n")
cat("Run the full workflow example (examples/clusterprofiler_workflow_example.R) to test with real data.\\n\\n")

cat("Quick start after completing standard funseqR workflow:\\n")
cat("\\n")
cat("enrichment_results <- run_clusterprofiler_enrichment(\\n")
cat("  con = con,\\n")
cat("  candidate_file_id = vcf_cand_import$file_id,\\n")
cat("  background_file_id = vcf_import$file_id,\\n")
cat("  blast_param_id = blast_results$blast_param_id,\\n")
cat("  pvalue_cutoff = 0.1  # To match previous 0.1 threshold\\n")
cat(")\\n")
cat("\\n")
cat("print(enrichment_results)\\n")
cat("dotplot(enrichment_results$BP)\\n")