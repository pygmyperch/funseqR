#!/usr/bin/env Rscript

# Test ORA function step by step with simplified architecture (no blast_param_id)
library(funseqR)

cat("=== Testing Simplified ORA Function (One Database = One Analysis) ===\n")

# Connect to database
con <- connect_funseq_db('funseq_project.db')

# Step 1: Test GO extraction helper function directly
cat("\n1. Testing .extract_stored_candidate_go_annotations()...\n")
tryCatch({
  go_result <- funseqR:::.extract_stored_candidate_go_annotations(con, verbose = TRUE)
  cat("SUCCESS: GO extraction completed\n")
  cat("  - Foreground genes:", length(go_result$foreground$genes), "\n")
  cat("  - Background genes:", length(go_result$background$genes), "\n")
  cat("  - Total GO terms:", nrow(go_result$all_go_terms), "\n")
}, error = function(e) {
  cat("ERROR in GO extraction:", e$message, "\n")
  close_funseq_db(con)
  quit(status = 1)
})

# Step 2: Test KEGG extraction helper function directly  
cat("\n2. Testing .extract_stored_candidate_kegg_annotations()...\n")
tryCatch({
  kegg_result <- funseqR:::.extract_stored_candidate_kegg_annotations(con, verbose = TRUE)
  cat("SUCCESS: KEGG extraction completed\n")
  cat("  - Foreground genes:", length(kegg_result$foreground$genes), "\n")
  cat("  - Background genes:", length(kegg_result$background$genes), "\n")
  cat("  - Total pathways:", nrow(kegg_result$all_pathways), "\n")
}, error = function(e) {
  cat("ERROR in KEGG extraction:", e$message, "\n")
  close_funseq_db(con)
  quit(status = 1)
})

# Step 3: Test main GO extraction function
cat("\n3. Testing extract_go_terms_for_enrichment()...\n")
tryCatch({
  go_data <- extract_go_terms_for_enrichment(con, "stored", verbose = TRUE)
  cat("SUCCESS: Main GO extraction function completed\n")
  cat("  - Foreground genes:", length(go_data$foreground$genes), "\n")
  cat("  - Background genes:", length(go_data$background$genes), "\n")
}, error = function(e) {
  cat("ERROR in main GO extraction:", e$message, "\n")
  close_funseq_db(con)
  quit(status = 1)
})

# Step 4: Test main KEGG extraction function
cat("\n4. Testing extract_kegg_terms_for_enrichment()...\n")
tryCatch({
  kegg_data <- extract_kegg_terms_for_enrichment(con, "stored", verbose = TRUE)
  cat("SUCCESS: Main KEGG extraction function completed\n")
  cat("  - Foreground genes:", length(kegg_data$foreground$genes), "\n")
  cat("  - Background genes:", length(kegg_data$background$genes), "\n")
}, error = function(e) {
  cat("ERROR in main KEGG extraction:", e$message, "\n")
  close_funseq_db(con)
  quit(status = 1)
})

# Step 5: Test full ORA function (GO only, no storage)
cat("\n5. Testing ora() function (GO only, no storage)...\n")
tryCatch({
  ora_result <- ora(con, 
                   candidate_vcf_file = "stored",
                   annotation_type = "GO",
                   significance_threshold = 0.1,
                   store_results = FALSE,  # Don't store yet
                   create_plots = FALSE,   # Don't create plots yet
                   verbose = TRUE)
  cat("SUCCESS: ORA function completed without storage\n")
  cat("  - Status:", ora_result$status, "\n")
  if (!is.null(ora_result$enrichment_results) && !is.null(ora_result$enrichment_results$GO)) {
    go_enrich <- ora_result$enrichment_results$GO
    if (!is.null(go_enrich$BP)) cat("  - BP enrichments:", nrow(go_enrich$BP), "\n")
    if (!is.null(go_enrich$MF)) cat("  - MF enrichments:", nrow(go_enrich$MF), "\n")
    if (!is.null(go_enrich$CC)) cat("  - CC enrichments:", nrow(go_enrich$CC), "\n")
  }
}, error = function(e) {
  cat("ERROR in ORA function:", e$message, "\n")
  close_funseq_db(con)
  quit(status = 1)
})

# Step 6: Test full ORA function with storage
cat("\n6. Testing ora() function with database storage...\n")
tryCatch({
  ora_result_stored <- ora(con, 
                          candidate_vcf_file = "stored",
                          annotation_type = "GO",
                          significance_threshold = 0.1,
                          store_results = TRUE,   # Test storage
                          create_plots = FALSE,  # Skip plots
                          verbose = TRUE)
  cat("SUCCESS: ORA function completed WITH storage\n")
  cat("  - Status:", ora_result_stored$status, "\n")
  if (!is.null(ora_result_stored$analysis_ids)) {
    cat("  - Analysis IDs stored:", paste(unlist(ora_result_stored$analysis_ids), collapse = ", "), "\n")
  }
}, error = function(e) {
  cat("ERROR in ORA function with storage:", e$message, "\n")
  close_funseq_db(con)
  quit(status = 1)
})

cat("\n=== All tests completed successfully! ===\n")
cat("The simplified ORA architecture (one database = one analysis) is working!\n")

# Close connection
close_funseq_db(con)