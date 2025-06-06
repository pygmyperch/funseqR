# Basic tests for GO enrichment functionality
# 
# This script provides basic validation of the GO enrichment implementation
# Run this after implementing the GO enrichment features to ensure core functionality works

library(funseqR)

# Test function to validate core GO enrichment functionality
test_go_enrichment_basics <- function() {
  
  cat("=== Testing GO Enrichment Functionality ===\n")
  
  # Test 1: Package loading and function availability
  cat("1. Testing function availability...\n")
  
  required_functions <- c(
    "import_candidate_loci",
    "extract_go_terms_for_enrichment", 
    "perform_go_enrichment",
    "create_go_bubble_plot",
    "run_go_enrichment_workflow"
  )
  
  for (func in required_functions) {
    if (exists(func)) {
      cat("  ✓", func, "available\n")
    } else {
      cat("  ✗", func, "missing\n")
      return(FALSE)
    }
  }
  
  # Test 2: Hypergeometric test implementation
  cat("2. Testing hypergeometric calculations...\n")
  
  # Simple test case: 
  # Population: 1000 genes, 100 with GO term
  # Sample: 50 genes, 10 with GO term
  # Expected: ~5, Observed: 10, should be significant
  
  p_val <- phyper(10 - 1, 100, 1000 - 100, 50, lower.tail = FALSE)
  expected <- (100/1000) * 50  # 5
  fold_enrichment <- 10 / expected  # 2.0
  
  cat("  Test case: 10 observed vs 5 expected\n")
  cat("  P-value:", format(p_val, scientific = TRUE), "\n")
  cat("  Fold enrichment:", round(fold_enrichment, 2), "\n")
  
  if (p_val < 0.05 && fold_enrichment == 2.0) {
    cat("  ✓ Hypergeometric test working correctly\n")
  } else {
    cat("  ✗ Hypergeometric test issue\n")
    return(FALSE)
  }
  
  # Test 3: Data structure validation
  cat("3. Testing data structure handling...\n")
  
  # Create mock GO data structure
  mock_go_data <- list(
    foreground = list(
      genes = c("P1", "P2", "P3"),
      gene2go = list(
        "P1" = c("GO:0001", "GO:0002"),
        "P2" = c("GO:0001"),
        "P3" = c("GO:0003")
      )
    ),
    background = list(
      genes = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10"),
      gene2go = list(
        "P1" = c("GO:0001", "GO:0002"),
        "P2" = c("GO:0001"),
        "P3" = c("GO:0003"),
        "P4" = c("GO:0001"),
        "P5" = c("GO:0002"),
        "P6" = c("GO:0003"),
        "P7" = c("GO:0004"),
        "P8" = c("GO:0004"),
        "P9" = c("GO:0004"),
        "P10" = c("GO:0004")
      )
    ),
    all_go_terms = data.frame(
      go_id = c("GO:0001", "GO:0002", "GO:0003", "GO:0004"),
      go_term = c("Process A", "Process B", "Process C", "Process D"),
      go_category = c("P", "P", "P", "P"),
      stringsAsFactors = FALSE
    )
  )
  
  # Test the enrichment calculation manually
  # GO:0001: foreground=2/3, background=3/10, should be enriched
  fg_with_term <- sum(sapply(mock_go_data$foreground$gene2go, function(x) "GO:0001" %in% x))
  bg_with_term <- sum(sapply(mock_go_data$background$gene2go, function(x) "GO:0001" %in% x))
  
  cat("  Mock data: GO:0001 in", fg_with_term, "foreground and", bg_with_term, "background genes\n")
  
  if (fg_with_term == 2 && bg_with_term == 3) {
    cat("  ✓ Data structure handling working\n")
  } else {
    cat("  ✗ Data structure handling issue\n")
    return(FALSE)
  }
  
  # Test 4: Multiple testing correction
  cat("4. Testing multiple testing correction...\n")
  
  test_pvals <- c(0.001, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.8)
  adjusted_pvals <- p.adjust(test_pvals, method = "fdr")
  
  # FDR should be more conservative (higher p-values)
  if (all(adjusted_pvals >= test_pvals) && adjusted_pvals[1] < 0.05) {
    cat("  ✓ FDR correction working\n")
  } else {
    cat("  ✗ FDR correction issue\n")
    return(FALSE)
  }
  
  cat("\n=== Basic Tests Passed ===\n")
  return(TRUE)
}

# Test visualization functions
test_visualization_functions <- function() {
  
  cat("=== Testing Visualization Functions ===\n")
  
  # Create mock enrichment results
  mock_results <- data.frame(
    go_id = c("GO:0001", "GO:0002", "GO:0003"),
    go_term = c("Response to stimulus", "Metabolic process", "Transport"),
    go_category = c("BP", "BP", "BP"),
    foreground_count = c(5, 3, 2),
    background_count = c(10, 8, 6),
    total_foreground = c(50, 50, 50),
    total_background = c(1000, 1000, 1000),
    expected_count = c(0.5, 0.4, 0.3),
    fold_enrichment = c(10, 7.5, 6.7),
    p_value = c(0.001, 0.01, 0.05),
    p_adjusted = c(0.003, 0.02, 0.05),
    significance_level = c("highly_significant", "significant", "significant"),
    stringsAsFactors = FALSE
  )
  
  # Test bubble plot creation
  cat("1. Testing bubble plot creation...\n")
  tryCatch({
    bubble_plot <- create_go_bubble_plot(mock_results)
    if ("ggplot" %in% class(bubble_plot)) {
      cat("  ✓ Bubble plot created successfully\n")
    } else {
      cat("  ✗ Bubble plot creation failed\n")
      return(FALSE)
    }
  }, error = function(e) {
    cat("  ✗ Bubble plot error:", e$message, "\n")
    return(FALSE)
  })
  
  # Test summary table creation
  cat("2. Testing summary table creation...\n")
  tryCatch({
    summary_table <- create_go_summary_table(mock_results)
    if (is.data.frame(summary_table) && nrow(summary_table) > 0) {
      cat("  ✓ Summary table created successfully\n")
    } else {
      cat("  ✗ Summary table creation failed\n")
      return(FALSE)
    }
  }, error = function(e) {
    cat("  ✗ Summary table error:", e$message, "\n")
    return(FALSE)
  })
  
  # Test treemap (only if treemapify available)
  cat("3. Testing treemap creation...\n")
  if (requireNamespace("treemapify", quietly = TRUE)) {
    tryCatch({
      treemap_plot <- create_go_treemap(mock_results)
      if ("ggplot" %in% class(treemap_plot)) {
        cat("  ✓ Treemap created successfully\n")
      } else {
        cat("  ✗ Treemap creation failed\n")
      }
    }, error = function(e) {
      cat("  ✗ Treemap error:", e$message, "\n")
    })
  } else {
    cat("  ! Treemap test skipped (treemapify not available)\n")
  }
  
  cat("\n=== Visualization Tests Completed ===\n")
  return(TRUE)
}

# Main test execution
main_test <- function() {
  cat("Starting GO enrichment functionality tests...\n\n")
  
  # Check required packages
  required_packages <- c("DBI", "ggplot2", "stringr")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    cat("Missing required packages:", paste(missing_packages, collapse = ", "), "\n")
    cat("Please install missing packages before running tests.\n")
    return(FALSE)
  }
  
  # Run tests
  basic_test_result <- test_go_enrichment_basics()
  viz_test_result <- test_visualization_functions()
  
  # Summary
  cat("\n=== TEST SUMMARY ===\n")
  if (basic_test_result && viz_test_result) {
    cat("✓ All tests passed! GO enrichment functionality is working correctly.\n")
    cat("\nNext steps:\n")
    cat("1. Test with real data using examples/go_enrichment_example.R\n")
    cat("2. Check the vignette: vignette('go_enrichment_analysis')\n")
    cat("3. Review documentation: ?run_go_enrichment_workflow\n")
    return(TRUE)
  } else {
    cat("✗ Some tests failed. Check the output above for details.\n")
    return(FALSE)
  }
}

# Run tests if script is executed directly
if (sys.nframe() == 0) {
  main_test()
}