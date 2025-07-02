#!/usr/bin/env Rscript

# Test script to verify the enrichment extraction fix
# This script tests the updated .extract_enrichment_data function

library(DBI)
library(dplyr)

# Source the updated function (assuming you're in the funseqR directory)
source("R/process_annotations.R")

cat("=== Testing Enrichment Extraction Fix ===\n")

# Connect to database (adjust path as needed)
# con <- dbConnect(RSQLite::SQLite(), "your_database.db")  # Replace with actual path

# Test function (replace with your actual database connection)
test_enrichment_extraction <- function(con) {
  
  cat("\n1. Testing with analysis_ids = c(1, 2, 3) and threshold = 0.1\n")
  
  # Test the updated function
  enrichment_data <- .extract_enrichment_data(con, c(1, 2, 3), 0.1, verbose = TRUE)
  
  cat("\nResults from updated function:\n")
  cat("Number of rows:", nrow(enrichment_data), "\n")
  
  if (nrow(enrichment_data) > 0) {
    cat("Columns:", paste(names(enrichment_data), collapse = ", "), "\n")
    cat("Analysis IDs found:", paste(unique(enrichment_data$analysis_id), collapse = ", "), "\n")
    cat("FDR range:", min(enrichment_data$fdr), "to", max(enrichment_data$fdr), "\n")
    
    # Show first few results
    cat("\nFirst few results:\n")
    print(head(enrichment_data[, c("analysis_id", "term_id", "term_name", "fdr")], 3))
  }
  
  cat("\n2. Testing edge cases\n")
  
  # Test with non-existent analysis IDs
  cat("\nTesting with non-existent analysis IDs (999, 1000):\n")
  empty_result <- .extract_enrichment_data(con, c(999, 1000), 0.1, verbose = TRUE)
  cat("Result rows:", nrow(empty_result), "\n")
  
  # Test with very strict threshold
  cat("\nTesting with very strict threshold (0.001):\n")
  strict_result <- .extract_enrichment_data(con, c(1, 2, 3), 0.001, verbose = TRUE)
  cat("Result rows:", nrow(strict_result), "\n")
  
  # Test with very lenient threshold
  cat("\nTesting with very lenient threshold (1.0):\n")
  lenient_result <- .extract_enrichment_data(con, c(1, 2, 3), 1.0, verbose = TRUE)
  cat("Result rows:", nrow(lenient_result), "\n")
  
  cat("\n3. Comparing with original query approach\n")
  
  # Test the original query for comparison
  analysis_ids_str <- paste(c(1, 2, 3), collapse = ", ")
  original_query <- paste0("
    SELECT DISTINCT
      ora.analysis_id,
      ora.annotation_type,
      ora.term_type,
      res.term_id,
      res.term_name,
      res.p_value,
      res.p_adjusted as fdr,
      res.gene_ids,
      res.fold_enrichment
    FROM ora_analyses ora
    JOIN ora_results res ON ora.analysis_id = res.analysis_id
    WHERE ora.analysis_id IN (", analysis_ids_str, ")
      AND res.p_adjusted <= ?
    ORDER BY res.p_adjusted
  ")
  
  original_result <- DBI::dbGetQuery(con, original_query, list(0.1))
  cat("Original query result rows:", nrow(original_result), "\n")
  cat("Updated function result rows:", nrow(enrichment_data), "\n")
  
  if (nrow(original_result) != nrow(enrichment_data)) {
    cat("DIFFERENCE DETECTED! The fix changed the results.\n")
    cat("This suggests the original issue was likely a data type problem.\n")
  } else {
    cat("Results are the same. The issue might be elsewhere.\n")
  }
  
  return(list(
    updated_result = enrichment_data,
    original_result = original_result,
    strict_result = strict_result,
    lenient_result = lenient_result
  ))
}

# Integration test with compile_funseq_results
test_integration <- function(con) {
  
  cat("\n=== Integration Test ===\n")
  
  # Test the full enrichment stage workflow
  cat("Testing compile_funseq_results enrichment stage...\n")
  
  # You would need to have annotation data from a previous stage
  # This is a placeholder - replace with actual data structure
  mock_annotation_data <- data.frame(
    locus_id = c("1_chr1_100", "1_chr1_200", "1_chr2_300"),
    uniprot_accession = c("P12345", "Q67890", "R11111"),
    gene_name = c("GENE1", "GENE2", "GENE3"),
    stringsAsFactors = FALSE
  )
  
  tryCatch({
    result <- compile_funseq_results(
      con = con,
      stage = "enrichment",
      data = mock_annotation_data,
      analysis_ids = c(1, 2, 3),
      significance_threshold = 0.1,
      verbose = TRUE
    )
    
    cat("Integration test completed successfully!\n")
    cat("Result rows:", nrow(result), "\n")
    cat("Enrichment columns present:", any(grepl("enriched", names(result))), "\n")
    
  }, error = function(e) {
    cat("Integration test failed with error:", e$message, "\n")
  })
}

# Manual testing instructions
cat("\n=== MANUAL TESTING INSTRUCTIONS ===\n")
cat("1. Connect to your database:\n")
cat("   con <- dbConnect(RSQLite::SQLite(), 'path_to_your_database.db')\n\n")
cat("2. Run the test function:\n")
cat("   test_results <- test_enrichment_extraction(con)\n\n")
cat("3. Run the integration test:\n")
cat("   test_integration(con)\n\n")
cat("4. Check the debug output to understand what changed\n\n")

cat("=== EXPECTED OUTCOMES ===\n")
cat("If the fix worked:\n")
cat("- You should now see enriched terms where there were none before\n")
cat("- The debug output will show the actual p_adjusted values and counts\n")
cat("- The difference between original and updated queries will be evident\n\n")

cat("If the issue persists:\n")
cat("- Check the debug output for clues about data types\n")
cat("- Verify that the analysis_ids exist in the database\n")
cat("- Look at the actual p_adjusted values to ensure they're within threshold\n")

cat("\n=== END TEST SCRIPT ===\n")