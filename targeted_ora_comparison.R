#!/usr/bin/env Rscript

# Targeted ORA Comparison: Focus on enrichment methodology differences
# Comparing cp_enrich_funseq_project.db vs simple_funseq_project.db

library(funseqR)
library(DBI)

cat("=== TARGETED ORA ENRICHMENT METHODOLOGY COMPARISON ===\n")
cat("Investigating differences in enrichment analysis methods\n\n")

# Connect to both databases
con_cp <- connect_funseq_db('cp_enrich_funseq_project.db')
con_simple <- connect_funseq_db('simple_funseq_project.db')

# Function to analyze ORA methodology
analyze_ora_methodology <- function(con, db_name) {
  cat("Analyzing", db_name, "methodology:\n")
  cat(paste(rep("-", 40), collapse = ""), "\n")
  
  # Check ora_analyses table for parameters used
  if ("ora_analyses" %in% DBI::dbListTables(con)) {
    analyses <- DBI::dbGetQuery(con, "
      SELECT analysis_id, annotation_type, term_type, 
             total_foreground_genes, total_background_genes,
             analysis_parameters, enrichment_method, analysis_date
      FROM ora_analyses 
      ORDER BY analysis_id
    ")
    
    cat("ORA Analyses Found:\n")
    if (nrow(analyses) > 0) {
      for (i in 1:nrow(analyses)) {
        row <- analyses[i, ]
        cat("  Analysis", row$analysis_id, ":\n")
        cat("    - Type:", row$annotation_type, "/", row$term_type, "\n")
        cat("    - Method:", row$enrichment_method, "\n")
        cat("    - Foreground genes:", row$total_foreground_genes, "\n")
        cat("    - Background genes:", row$total_background_genes, "\n")
        cat("    - Parameters:", row$analysis_parameters, "\n")
        cat("    - Date:", row$analysis_date, "\n")
      }
    } else {
      cat("  No analyses found\n")
    }
    cat("\n")
    return(analyses)
  } else {
    cat("  No ora_analyses table found\n\n")
    return(data.frame())
  }
}

# Function to compare specific enrichment results
compare_enrichment_results <- function(con_cp, con_simple) {
  cat("Comparing Enrichment Results:\n")
  cat(paste(rep("=", 40), collapse = ""), "\n")
  
  # Get results from both databases for comparison
  if ("ora_results" %in% DBI::dbListTables(con_cp) && 
      "ora_results" %in% DBI::dbListTables(con_simple)) {
    
    # Get GO:0036342 (post-anal tail morphogenesis) from both - this appears in both result sets
    cp_results <- DBI::dbGetQuery(con_cp, "
      SELECT analysis_id, term_id, term_name, 
             foreground_count, background_count, total_foreground, total_background,
             expected_count, fold_enrichment, p_value, p_adjusted,
             gene_ratio, bg_ratio, gene_ids
      FROM ora_results 
      WHERE term_id = 'GO:0036342'
      ORDER BY analysis_id
    ")
    
    simple_results <- DBI::dbGetQuery(con_simple, "
      SELECT analysis_id, term_id, term_name,
             foreground_count, background_count, total_foreground, total_background, 
             expected_count, fold_enrichment, p_value, p_adjusted,
             gene_ratio, bg_ratio, gene_ids
      FROM ora_results 
      WHERE term_id = 'GO:0036342'
      ORDER BY analysis_id
    ")
    
    cat("GO:0036342 (post-anal tail morphogenesis) comparison:\n")
    
    if (nrow(cp_results) > 0) {
      cat("  cp_enrich database:\n")
      for (i in 1:nrow(cp_results)) {
        row <- cp_results[i, ]
        cat("    - foreground_count:", row$foreground_count, "\n")
        cat("    - background_count:", row$background_count, "\n") 
        cat("    - total_foreground:", row$total_foreground, "\n")
        cat("    - total_background:", row$total_background, "\n")
        cat("    - p_value:", row$p_value, "\n")
        cat("    - p_adjusted:", row$p_adjusted, "\n")
        cat("    - fold_enrichment:", row$fold_enrichment, "\n")
        cat("    - gene_ratio:", row$gene_ratio, "\n")
        cat("    - bg_ratio:", row$bg_ratio, "\n")
        cat("    - gene_ids:", row$gene_ids, "\n")
      }
    } else {
      cat("  cp_enrich: GO:0036342 not found\n")
    }
    
    if (nrow(simple_results) > 0) {
      cat("  simple database:\n")
      for (i in 1:nrow(simple_results)) {
        row <- simple_results[i, ]
        cat("    - foreground_count:", row$foreground_count, "\n")
        cat("    - background_count:", row$background_count, "\n")
        cat("    - total_foreground:", row$total_foreground, "\n") 
        cat("    - total_background:", row$total_background, "\n")
        cat("    - p_value:", row$p_value, "\n")
        cat("    - p_adjusted:", row$p_adjusted, "\n")
        cat("    - fold_enrichment:", row$fold_enrichment, "\n")
        cat("    - gene_ratio:", row$gene_ratio, "\n")
        cat("    - bg_ratio:", row$bg_ratio, "\n")
        cat("    - gene_ids:", row$gene_ids, "\n")
      }
    } else {
      cat("  simple: GO:0036342 not found\n")
    }
    
    cat("\n")
    
    # Also check P:hemopoiesis (GO:0030097) which appears in cp_enrich results
    cat("GO:0030097 (hemopoiesis) comparison:\n")
    
    cp_hemo <- DBI::dbGetQuery(con_cp, "
      SELECT analysis_id, term_id, term_name,
             foreground_count, background_count, total_foreground, total_background,
             p_value, p_adjusted, fold_enrichment, gene_ids
      FROM ora_results 
      WHERE term_id = 'GO:0030097'
    ")
    
    simple_hemo <- DBI::dbGetQuery(con_simple, "
      SELECT analysis_id, term_id, term_name,
             foreground_count, background_count, total_foreground, total_background,
             p_value, p_adjusted, fold_enrichment, gene_ids  
      FROM ora_results 
      WHERE term_id = 'GO:0030097'
    ")
    
    if (nrow(cp_hemo) > 0) {
      cat("  cp_enrich found GO:0030097:\n")
      row <- cp_hemo[1, ]
      cat("    - p_adjusted:", row$p_adjusted, " (significant in cp_enrich)\n")
      cat("    - foreground_count:", row$foreground_count, "\n")
      cat("    - gene_ids:", row$gene_ids, "\n")
    } else {
      cat("  cp_enrich: GO:0030097 not found\n")
    }
    
    if (nrow(simple_hemo) > 0) {
      cat("  simple found GO:0030097:\n") 
      row <- simple_hemo[1, ]
      cat("    - p_adjusted:", row$p_adjusted, "\n")
      cat("    - foreground_count:", row$foreground_count, "\n")
      cat("    - gene_ids:", row$gene_ids, "\n")
    } else {
      cat("  simple: GO:0030097 not found or not significant\n")
    }
  }
}

# Function to extract and compare gene sets
compare_gene_sets <- function(con_cp, con_simple) {
  cat("\nComparing Gene Set Extraction:\n")
  cat(paste(rep("=", 40), collapse = ""), "\n")
  
  tryCatch({
    # Extract GO data from both databases using current method
    cat("Extracting GO data using current implementation...\n")
    
    go_data_cp <- extract_go_terms_for_enrichment(con_cp, "stored", verbose = FALSE)
    go_data_simple <- extract_go_terms_for_enrichment(con_simple, "stored", verbose = FALSE)
    
    cat("Gene set sizes:\n")
    cat("  cp_enrich - foreground:", length(go_data_cp$foreground$genes), 
        "background:", length(go_data_cp$background$genes), "\n")
    cat("  simple - foreground:", length(go_data_simple$foreground$genes),
        "background:", length(go_data_simple$background$genes), "\n")
    
    # Check if gene sets are identical
    fg_same <- setequal(go_data_cp$foreground$genes, go_data_simple$foreground$genes)
    bg_same <- setequal(go_data_cp$background$genes, go_data_simple$background$genes)
    
    cat("Gene sets identical:\n")
    cat("  Foreground:", fg_same, "\n")
    cat("  Background:", bg_same, "\n")
    
    if (!fg_same) {
      fg_diff <- setdiff(go_data_cp$foreground$genes, go_data_simple$foreground$genes)
      if (length(fg_diff) > 0) {
        cat("  cp_enrich only foreground:", paste(head(fg_diff, 3), collapse = ", "), "\n")
      }
      fg_diff2 <- setdiff(go_data_simple$foreground$genes, go_data_cp$foreground$genes)
      if (length(fg_diff2) > 0) {
        cat("  simple only foreground:", paste(head(fg_diff2, 3), collapse = ", "), "\n")
      }
    }
    
    return(list(cp = go_data_cp, simple = go_data_simple))
    
  }, error = function(e) {
    cat("Error extracting gene sets:", e$message, "\n")
    return(NULL)
  })
}

# Run the analysis
cat("1. METHODOLOGY ANALYSIS\n")
cat("========================\n")

analyses_cp <- analyze_ora_methodology(con_cp, "cp_enrich")
analyses_simple <- analyze_ora_methodology(con_simple, "simple")

cat("2. ENRICHMENT RESULTS COMPARISON\n")
cat("================================\n")

compare_enrichment_results(con_cp, con_simple)

cat("3. GENE SET COMPARISON\n")
cat("======================\n")

gene_sets <- compare_gene_sets(con_cp, con_simple)

cat("4. SUMMARY AND CONCLUSIONS\n")
cat("==========================\n")

# Compare key parameters
if (nrow(analyses_cp) > 0 && nrow(analyses_simple) > 0) {
  cat("Key Differences Found:\n")
  
  # Check methods
  cp_methods <- unique(analyses_cp$enrichment_method)
  simple_methods <- unique(analyses_simple$enrichment_method)
  
  cat("• Enrichment methods:\n")
  cat("  - cp_enrich:", paste(cp_methods, collapse = ", "), "\n")
  cat("  - simple:", paste(simple_methods, collapse = ", "), "\n")
  
  if (!identical(cp_methods, simple_methods)) {
    cat("  ⚠️  DIFFERENT ENRICHMENT METHODS! This explains the different results.\n")
  }
  
  # Check gene set sizes
  cp_fg <- unique(analyses_cp$total_foreground_genes)
  simple_fg <- unique(analyses_simple$total_foreground_genes)
  cp_bg <- unique(analyses_cp$total_background_genes)
  simple_bg <- unique(analyses_simple$total_background_genes)
  
  cat("• Gene set sizes:\n")
  cat("  - cp_enrich foreground:", paste(cp_fg, collapse = ", "), "\n")
  cat("  - simple foreground:", paste(simple_fg, collapse = ", "), "\n")
  cat("  - cp_enrich background:", paste(cp_bg, collapse = ", "), "\n")
  cat("  - simple background:", paste(simple_bg, collapse = ", "), "\n")
  
  if (!identical(cp_fg, simple_fg) || !identical(cp_bg, simple_bg)) {
    cat("  ⚠️  DIFFERENT GENE SET SIZES! This could explain different results.\n")
  }
}

cat("\n=== ANALYSIS COMPLETE ===\n")

# Close connections
close_funseq_db(con_cp)
close_funseq_db(con_simple)