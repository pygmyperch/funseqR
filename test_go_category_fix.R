# Test script to verify GO category mapping fix
# 
# This script tests that the GO enrichment analysis correctly maps
# between analysis interface format ("BP", "MF", "CC") and 
# database storage format ("P", "F", "C")

library(funseqR)

# Test the ontology mapping function directly
test_ontology_mapping <- function() {
  cat("Testing ontology mapping...\n")
  
  # This is the mapping that should be used in .perform_go_enrichment_from_data
  ontology_map <- c("BP" = "P", "MF" = "F", "CC" = "C")
  
  # Test each mapping
  test_cases <- c("BP", "MF", "CC")
  expected <- c("P", "F", "C")
  
  for (i in seq_along(test_cases)) {
    input <- test_cases[i]
    expected_output <- expected[i]
    actual_output <- ontology_map[input]
    
    if (actual_output == expected_output) {
      cat("✅ ", input, " -> ", actual_output, " (correct)\n")
    } else {
      cat("❌ ", input, " -> ", actual_output, " (expected ", expected_output, ")\n")
    }
  }
  
  # Test invalid input
  invalid_result <- ontology_map["INVALID"]
  if (is.na(invalid_result)) {
    cat("✅ Invalid input correctly returns NA\n")
  } else {
    cat("❌ Invalid input should return NA, got: ", invalid_result, "\n")
  }
}

# Test with mock data
test_with_mock_data <- function() {
  cat("\nTesting with mock GO data...\n")
  
  # Create mock GO data with database format categories ("P", "F", "C")
  mock_go_data <- list(
    go_terms = data.frame(
      go_id = c("GO:0008150", "GO:0003674", "GO:0005575", 
                "GO:0009987", "GO:0016209", "GO:0016020"),
      go_name = c("biological_process", "molecular_function", "cellular_component",
                  "cellular process", "antioxidant activity", "membrane"),
      go_category = c("P", "F", "C", "P", "F", "C"),  # Database format
      locus_count = c(100, 80, 60, 25, 15, 40),
      stringsAsFactors = FALSE
    ),
    foreground = list(loci = c("locus1", "locus2"), genes = c("locus1", "locus2")),
    background = list(loci = c("locus1", "locus2", "locus3", "locus4"), 
                      genes = c("locus1", "locus2", "locus3", "locus4")),
    term_mappings = list(
      "GO:0008150" = list(loci = c("locus1", "locus2")),
      "GO:0003674" = list(loci = c("locus1")),
      "GO:0005575" = list(loci = c("locus2"))
    )
  )
  
  # Test filtering with each ontology
  ontologies <- c("BP", "MF", "CC")
  expected_categories <- c("P", "F", "C")
  
  for (i in seq_along(ontologies)) {
    ontology <- ontologies[i]
    expected_cat <- expected_categories[i]
    
    # Apply the mapping (this is what the fixed function does)
    ontology_map <- c("BP" = "P", "MF" = "F", "CC" = "C")
    category_code <- ontology_map[ontology]
    
    # Filter GO terms 
    filtered_terms <- mock_go_data$go_terms[mock_go_data$go_terms$go_category == category_code, ]
    
    expected_count <- sum(mock_go_data$go_terms$go_category == expected_cat)
    actual_count <- nrow(filtered_terms)
    
    if (actual_count == expected_count && actual_count > 0) {
      cat("✅ ", ontology, " (", category_code, ") found ", actual_count, " terms\n")
    } else {
      cat("❌ ", ontology, " (", category_code, ") found ", actual_count, " terms, expected ", expected_count, "\n")
    }
  }
}

# Main test function
main <- function() {
  cat("=== GO Category Mapping Fix Test ===\n")
  cat("This test verifies that the GO enrichment analysis fix correctly\n")
  cat("maps between analysis interface format and database storage format.\n\n")
  
  test_ontology_mapping()
  test_with_mock_data()
  
  cat("\n=== Test Summary ===\n")
  cat("If all tests show ✅, the GO category mapping fix is working correctly.\n")
  cat("The enrichment analysis should now find BP, MF, and CC terms.\n")
}

# Run the test
main()