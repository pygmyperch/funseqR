# Manual GO Enrichment Test
# This script manually performs Fisher's exact test to verify expected results

manual_go_enrichment_test <- function(annotations, candidate_loci, ontology = "BP", verbose = TRUE) {
  
  cat("=== Manual GO Enrichment Test ===\n")
  cat("Ontology:", ontology, "\n\n")
  
  # Step 1: Get candidate locus IDs 
  bed_coords <- candidate_loci
  chrom_col <- if ("chrom" %in% colnames(bed_coords)) "chrom" else "chromosome"
  
  candidate_df <- data.frame(
    chromosome = bed_coords[[chrom_col]],
    position = bed_coords$start + 1,
    stringsAsFactors = FALSE
  )
  
  candidate_patterns <- with(candidate_df, paste0("_", chromosome, "_", position, "$"))
  
  candidate_locus_ids <- annotations$locus_id[
    sapply(annotations$locus_id, function(locus_id) {
      any(sapply(candidate_patterns, function(pattern) {
        grepl(pattern, locus_id)
      }))
    })
  ]
  
  if (verbose) {
    cat("Step 1 - Candidate identification:\n")
    cat("  Candidate loci found:", length(candidate_locus_ids), "\n")
  }
  
  # Step 2: Filter annotations with GO terms
  go_annotations <- annotations[!is.na(annotations$go_terms) & annotations$go_terms != "", ]
  
  if (verbose) {
    cat("  Background with GO:", nrow(go_annotations), "\n")
  }
  
  # Step 3: Get candidates with GO annotations
  candidate_with_go <- intersect(candidate_locus_ids, go_annotations$locus_id)
  
  if (verbose) {
    cat("  Candidates with GO:", length(candidate_with_go), "\n\n")
  }
  
  if (length(candidate_with_go) == 0) {
    cat("ERROR: No candidates have GO annotations - cannot perform enrichment\n")
    return(NULL)
  }
  
  # Step 4: Extract GO terms and build term map
  if (verbose) cat("Step 2 - Building GO term map:\n")
  
  go_term_map <- list()
  
  for (i in 1:nrow(go_annotations)) {
    locus_id <- go_annotations$locus_id[i]
    
    if (!is.na(go_annotations$go_terms[i]) && go_annotations$go_terms[i] != "") {
      go_terms <- unlist(strsplit(go_annotations$go_terms[i], ";"))
      go_categories <- unlist(strsplit(go_annotations$go_categories[i], ";"))
      
      for (j in seq_along(go_terms)) {
        go_id <- go_terms[j]
        go_category <- if (j <= length(go_categories)) go_categories[j] else "Unknown"
        
        if (!is.na(go_id) && go_id != "") {
          if (!go_id %in% names(go_term_map)) {
            go_term_map[[go_id]] <- list(
              go_id = go_id,
              go_category = go_category,
              loci = character(0)
            )
          }
          go_term_map[[go_id]]$loci <- c(go_term_map[[go_id]]$loci, locus_id)
        }
      }
    }
  }
  
  if (verbose) {
    cat("  Total GO terms found:", length(go_term_map), "\n")
  }
  
  # Step 5: Filter for specific ontology
  ontology_map <- c("BP" = "P", "MF" = "F", "CC" = "C")
  target_category <- ontology_map[ontology]
  
  ontology_terms <- go_term_map[sapply(go_term_map, function(x) x$go_category == target_category)]
  
  if (verbose) {
    cat("  Terms in", ontology, "ontology:", length(ontology_terms), "\n\n")
  }
  
  if (length(ontology_terms) == 0) {
    cat("ERROR: No", ontology, "terms found\n")
    return(NULL)
  }
  
  # Step 6: Perform Fisher's exact test for each term
  if (verbose) cat("Step 3 - Testing for enrichment:\n")
  
  results <- data.frame()
  tests_performed <- 0
  
  for (i in 1:min(length(ontology_terms), 20)) {  # Test first 20 terms
    term_info <- ontology_terms[[i]]
    go_id <- term_info$go_id
    term_loci <- unique(term_info$loci)
    
    # Skip terms with too few or too many annotations
    if (length(term_loci) < 3 || length(term_loci) > 500) next
    
    # Calculate Fisher's exact test components
    candidate_with_term <- length(intersect(candidate_locus_ids, term_loci))
    candidate_without_term <- length(candidate_locus_ids) - candidate_with_term
    background_with_term <- length(intersect(go_annotations$locus_id, term_loci)) - candidate_with_term
    background_without_term <- length(go_annotations$locus_id) - length(intersect(go_annotations$locus_id, term_loci)) - candidate_without_term
    
    # Only test if we have candidates with this term
    if (candidate_with_term > 0) {
      tests_performed <- tests_performed + 1
      
      # Create contingency table
      cont_table <- matrix(c(candidate_with_term, candidate_without_term,
                            background_with_term, background_without_term), 
                          nrow = 2, byrow = TRUE)
      
      # Perform Fisher's exact test
      test_result <- fisher.test(cont_table, alternative = "greater")
      
      # Calculate enrichment metrics
      total_background <- length(go_annotations$locus_id)
      total_candidates <- length(candidate_locus_ids)
      
      fold_enrichment <- (candidate_with_term / total_candidates) / 
                        (length(intersect(go_annotations$locus_id, term_loci)) / total_background)
      
      result_row <- data.frame(
        go_id = go_id,
        candidate_count = candidate_with_term,
        background_count = length(intersect(go_annotations$locus_id, term_loci)),
        total_candidates = total_candidates,
        total_background = total_background,
        fold_enrichment = fold_enrichment,
        p_value = test_result$p.value,
        stringsAsFactors = FALSE
      )
      
      results <- rbind(results, result_row)
      
      if (verbose && test_result$p.value < 0.2) {  # Show promising results
        cat(sprintf("  %s: %d/%d candidates, %d/%d background, p=%.4f, fold=%.2f\n", 
                   go_id, candidate_with_term, total_candidates, 
                   length(intersect(go_annotations$locus_id, term_loci)), total_background,
                   test_result$p.value, fold_enrichment))
      }
    }
  }
  
  if (nrow(results) > 0) {
    # Adjust p-values
    results$p_adjusted <- p.adjust(results$p_value, method = "fdr")
    
    # Sort by p-value
    results <- results[order(results$p_value), ]
    
    if (verbose) {
      cat("\nStep 4 - Results summary:\n")
      cat("  Tests performed:", tests_performed, "\n")
      cat("  Significant at p < 0.05:", sum(results$p_value < 0.05), "\n") 
      cat("  Significant at FDR < 0.1:", sum(results$p_adjusted < 0.1), "\n")
      
      if (sum(results$p_adjusted < 0.1) > 0) {
        cat("\nSignificant results (FDR < 0.1):\n")
        sig_results <- results[results$p_adjusted < 0.1, ]
        for (i in 1:nrow(sig_results)) {
          cat(sprintf("  %s: p=%.4f, FDR=%.4f, fold=%.2f\n",
                     sig_results$go_id[i], sig_results$p_value[i], 
                     sig_results$p_adjusted[i], sig_results$fold_enrichment[i]))
        }
      }
    }
  } else {
    if (verbose) {
      cat("  No terms had candidates - cannot perform enrichment test\n")
    }
  }
  
  cat("\n=== Manual Test Complete ===\n")
  return(results)
}

# Usage:
# manual_results <- manual_go_enrichment_test(annotations, candidates$bed_file, "BP")