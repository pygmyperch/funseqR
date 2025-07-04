# Debug script to trace clusterProfiler values through the pipeline
# This will help identify exactly where p_adjusted gets corrupted

debug_clusterprofiler_conversion <- function() {
  
  cat("=== Debugging clusterProfiler p_adjusted Conversion ===\n")
  
  # We need to monkey-patch the conversion function to add debug output
  # First, let's save the original function
  original_convert <- funseqR:::.convert_clusterprofiler_to_funseqr
  
  # Create a debug version
  debug_convert <- function(clusterprofiler_result, ontology, significance_threshold, verbose) {
    
    cat("🔍 DEBUG: Entering .convert_clusterprofiler_to_funseqr\n")
    
    if (is.null(clusterprofiler_result) || nrow(clusterprofiler_result@result) == 0) {
      if (verbose) message("  - No terms found by clusterProfiler")
      return(data.frame())
    }
    
    cp_df <- clusterprofiler_result@result
    
    cat("🔍 DEBUG: clusterProfiler raw results:\n")
    cat("  - Number of terms:", nrow(cp_df), "\n")
    if (nrow(cp_df) > 0) {
      cat("  - p.adjust column class:", class(cp_df$p.adjust), "\n")
      cat("  - p.adjust values (first 5):", head(cp_df$p.adjust, 5), "\n")
      cat("  - p.adjust range:", range(cp_df$p.adjust, na.rm = TRUE), "\n")
    }
    
    # Extract numeric values from ratios for calculations
    fg_count <- cp_df$Count
    bg_count <- as.numeric(sub("/.*", "", cp_df$BgRatio))
    total_fg <- as.numeric(sub(".*/", "", cp_df$GeneRatio))
    total_bg <- as.numeric(sub(".*/", "", cp_df$BgRatio))
    
    # Calculate expected count and fold enrichment correctly
    expected_count <- (bg_count / total_bg) * total_fg
    fold_enrichment <- ifelse(expected_count > 0, fg_count / expected_count, Inf)
    
    cat("🔍 DEBUG: After calculations:\n")
    cat("  - fg_count range:", range(fg_count, na.rm = TRUE), "\n")
    cat("  - expected_count range:", range(expected_count, na.rm = TRUE), "\n")
    
    # Convert p_adjusted explicitly
    p_adjusted_numeric <- as.numeric(cp_df$p.adjust)
    cat("🔍 DEBUG: p_adjusted conversion:\n")
    cat("  - Original p.adjust:", head(cp_df$p.adjust, 3), "\n")
    cat("  - as.numeric() result:", head(p_adjusted_numeric, 3), "\n")
    cat("  - Any NAs introduced?", any(is.na(p_adjusted_numeric)), "\n")
    
    # Convert to funseqR format with additional clusterProfiler fields
    funseqr_results <- data.frame(
      go_id = cp_df$ID,
      go_term = cp_df$Description,
      go_category = ontology,
      foreground_count = fg_count,
      background_count = bg_count,
      total_foreground = total_fg,
      total_background = total_bg,
      expected_count = expected_count,
      fold_enrichment = fold_enrichment,
      p_value = as.numeric(cp_df$pvalue),
      p_adjusted = p_adjusted_numeric,
      significance_level = ifelse(cp_df$p.adjust < (significance_threshold / 5), "highly_significant",
                                 ifelse(cp_df$p.adjust < significance_threshold, "significant",
                                       ifelse(cp_df$p.adjust < (significance_threshold * 2), "trending", "not_significant"))),
      gene_ratio = cp_df$GeneRatio,
      bg_ratio = cp_df$BgRatio,
      qvalue = as.numeric(cp_df$qvalue),
      gene_ids = cp_df$geneID,
      stringsAsFactors = FALSE
    )
    
    cat("🔍 DEBUG: Final funseqr_results:\n")
    cat("  - p_adjusted column class:", class(funseqr_results$p_adjusted), "\n")
    cat("  - p_adjusted values (first 3):", head(funseqr_results$p_adjusted, 3), "\n")
    cat("  - p_adjusted range:", range(funseqr_results$p_adjusted, na.rm = TRUE), "\n")
    
    # Check significance
    sig_count <- sum(funseqr_results$p_adjusted < significance_threshold, na.rm = TRUE)
    cat("  - Significant terms (< ", significance_threshold, "):", sig_count, "\n")
    
    if (verbose) {
      message("  - clusterProfiler analysis complete: ", nrow(funseqr_results), " terms tested, ", sig_count, " significantly enriched (FDR < ", significance_threshold, ")")
    }
    
    cat("🔍 DEBUG: Exiting .convert_clusterprofiler_to_funseqr\n")
    return(funseqr_results)
  }
  
  # Replace the function in the namespace
  environment(debug_convert) <- environment(original_convert)
  assignInNamespace(".convert_clusterprofiler_to_funseqr", debug_convert, "funseqR")
  
  cat("✅ Debug version installed. Re-run your ORA analysis to see debug output.\n")
  cat("   Use this command to restore original function:\n")
  cat("   assignInNamespace('.convert_clusterprofiler_to_funseqr', original_convert, 'funseqR')\n")
  
  # Return the original function so user can restore it
  return(original_convert)
}

# Run the debug setup
original_function <- debug_clusterprofiler_conversion()

cat("\n🔧 Debug setup complete! Now re-run a simple ORA test:\n")
cat("   delete_ora_results(con, analysis_id = c(1, 2, 3))\n") 
cat("   ORA_results <- run_ORA(con, vcf_file_cand, annotation_type = 'GO', significance_threshold = 0.1, blast_param_id = 1)\n")