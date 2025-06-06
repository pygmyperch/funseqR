#' GO Enrichment Summary Template Functions
#'
#' Functions to generate formatted text summaries of GO enrichment results
#' for reports and standalone summary files.

#' Generate GO enrichment text summary
#'
#' @param ontology_summary Data frame with ontology testing overview
#' @param enriched_terms Data frame with enriched GO terms
#' @param gene_count Integer. Number of genes contributing to enrichment
#' @param output_file Character. Optional file path to write summary
#' @param format Character. Output format: "text" or "markdown". Default "text"
#'
#' @return Character vector containing formatted summary text
#'
#' @export
generate_go_summary_text <- function(ontology_summary, enriched_terms, gene_count, 
                                    output_file = NULL, format = "text") {
  
  # Initialize summary lines
  summary_lines <- c()
  
  # Header
  if (format == "markdown") {
    summary_lines <- c(summary_lines, "## GO Enrichment Analysis Summary", "")
  } else {
    summary_lines <- c(summary_lines, "GO Enrichment Analysis Summary", 
                      paste(rep("=", 35), collapse = ""), "")
  }
  
  # Testing Overview
  if (format == "markdown") {
    summary_lines <- c(summary_lines, "### Testing Overview:")
  } else {
    summary_lines <- c(summary_lines, "Testing Overview:")
  }
  summary_lines <- c(summary_lines, "")
  
  for (i in 1:nrow(ontology_summary)) {
    ont <- ontology_summary[i, ]
    ont_name <- switch(ont$ontology,
      "BP" = "Biological Process",
      "MF" = "Molecular Function", 
      "CC" = "Cellular Component"
    )
    
    line <- sprintf("  - %s: %d terms tested → %d enriched (FDR < 0.1)",
                   ont_name, ont$total_tested, ont$trending)
    summary_lines <- c(summary_lines, line)
  }
  
  summary_lines <- c(summary_lines, "")
  
  # Enriched Terms Section
  if (nrow(enriched_terms) > 0) {
    total_enriched <- nrow(enriched_terms)
    
    if (format == "markdown") {
      summary_lines <- c(summary_lines, sprintf("### Unique Enriched GO Terms (%d total):", total_enriched))
    } else {
      summary_lines <- c(summary_lines, sprintf("Unique Enriched GO Terms (%d total):", total_enriched))
    }
    summary_lines <- c(summary_lines, "")
    
    # Group by significance level
    highly_sig <- enriched_terms[enriched_terms$p_adjusted < 0.01, ]
    significant <- enriched_terms[enriched_terms$p_adjusted >= 0.01 & 
                                enriched_terms$p_adjusted < 0.05, ]
    trending <- enriched_terms[enriched_terms$p_adjusted >= 0.05 & 
                             enriched_terms$p_adjusted < 0.1, ]
    
    counter <- 1
    
    # Highly significant
    if (nrow(highly_sig) > 0) {
      if (format == "markdown") {
        summary_lines <- c(summary_lines, "**Highly Significant (FDR < 0.01):**")
      } else {
        summary_lines <- c(summary_lines, "Highly Significant (FDR < 0.01):")
      }
      
      for (i in 1:nrow(highly_sig)) {
        term <- highly_sig[i, ]
        line <- sprintf("%d. %s - %d genes (%.1f× enriched)", 
                       counter, term$go_term, term$foreground_count, term$fold_enrichment)
        summary_lines <- c(summary_lines, line)
        counter <- counter + 1
      }
      summary_lines <- c(summary_lines, "")
    }
    
    # Significant  
    if (nrow(significant) > 0) {
      if (format == "markdown") {
        summary_lines <- c(summary_lines, "**Significant (FDR < 0.05):**")
      } else {
        summary_lines <- c(summary_lines, "Significant (FDR < 0.05):")
      }
      
      for (i in 1:nrow(significant)) {
        term <- significant[i, ]
        line <- sprintf("%d. %s - %d genes (%.1f× enriched)",
                       counter, term$go_term, term$foreground_count, term$fold_enrichment)
        summary_lines <- c(summary_lines, line)
        counter <- counter + 1
      }
      summary_lines <- c(summary_lines, "")
    }
    
    # Trending
    if (nrow(trending) > 0) {
      if (format == "markdown") {
        summary_lines <- c(summary_lines, "**Trending (FDR < 0.1):**")
      } else {
        summary_lines <- c(summary_lines, "Trending (FDR < 0.1):")
      }
      
      for (i in 1:nrow(trending)) {
        term <- trending[i, ]
        line <- sprintf("%d. %s - %d genes (%.1f× enriched)",
                       counter, term$go_term, term$foreground_count, term$fold_enrichment)
        summary_lines <- c(summary_lines, line)
        counter <- counter + 1
      }
      summary_lines <- c(summary_lines, "")
    }
    
  } else {
    summary_lines <- c(summary_lines, "No significantly enriched GO terms found (FDR < 0.1)", "")
  }
  
  # Key Insights
  if (nrow(enriched_terms) > 0) {
    if (format == "markdown") {
      summary_lines <- c(summary_lines, "### Key Insights:")
    } else {
      summary_lines <- c(summary_lines, "Key Insights:")
    }
    summary_lines <- c(summary_lines, "")
    
    # Gene count insight
    summary_lines <- c(summary_lines, sprintf("  - %d unique genes contribute to the enrichment signal", gene_count))
    
    # Ontology-specific insights
    bp_count <- sum(ontology_summary$ontology == "BP" & ontology_summary$trending > 0)
    mf_count <- sum(ontology_summary$ontology == "MF" & ontology_summary$trending > 0) 
    cc_count <- sum(ontology_summary$ontology == "CC" & ontology_summary$trending > 0)
    
    if (bp_count == 0) {
      summary_lines <- c(summary_lines, "  - No Biological Process terms are significantly enriched")
    }
    
    # Analyze enriched terms for functional insights
    if (nrow(enriched_terms) > 0) {
      
      # Look for localization patterns
      localization_terms <- c("perikaryon", "extracellular region", "nucleus", "cytoplasm", 
                             "membrane", "organelle")
      has_localization <- any(sapply(localization_terms, function(x) 
        any(grepl(x, enriched_terms$go_term, ignore.case = TRUE))))
      
      if (has_localization) {
        loc_terms <- enriched_terms[grepl(paste(localization_terms, collapse = "|"), 
                                        enriched_terms$go_term, ignore.case = TRUE), ]
        if (nrow(loc_terms) > 0) {
          loc_names <- gsub("^[CFP]:", "", loc_terms$go_term)
          summary_lines <- c(summary_lines, 
            sprintf("  - Strong enrichment in cellular localization (%s)", 
                   paste(loc_names, collapse = ", ")))
        }
      }
      
      # Look for binding patterns
      binding_terms <- enriched_terms[grepl("binding", enriched_terms$go_term, ignore.case = TRUE), ]
      if (nrow(binding_terms) > 0) {
        binding_names <- gsub("^[CFP]:", "", binding_terms$go_term)
        summary_lines <- c(summary_lines,
          sprintf("  - Moderate enrichment in binding functions (%s)", 
                 paste(binding_names, collapse = ", ")))
      }
      
      # Look for metabolic patterns
      metabolic_terms <- enriched_terms[grepl("metabolic|synthesis|degradation|catalytic", 
                                             enriched_terms$go_term, ignore.case = TRUE), ]
      if (nrow(metabolic_terms) > 0) {
        metab_names <- gsub("^[CFP]:", "", metabolic_terms$go_term)
        summary_lines <- c(summary_lines,
          sprintf("  - Enrichment in metabolic processes (%s)", 
                 paste(metab_names, collapse = ", ")))
      }
    }
  }
  
  # Write to file if requested
  if (!is.null(output_file)) {
    writeLines(summary_lines, output_file)
    message("Summary written to: ", output_file)
  }
  
  return(summary_lines)
}

#' Generate GO enrichment summary for dynamic reports
#'
#' @param ontology_summary Data frame with ontology testing overview
#' @param enriched_terms Data frame with enriched GO terms  
#' @param gene_count Integer. Number of genes contributing to enrichment
#' @param analysis_date Character. Date of analysis
#'
#' @return List with summary sections for dynamic report integration
#'
#' @export
generate_go_report_content <- function(ontology_summary, enriched_terms, gene_count, 
                                      analysis_date = Sys.Date()) {
  
  # Generate markdown summary
  summary_text <- generate_go_summary_text(ontology_summary, enriched_terms, 
                                          gene_count, format = "markdown")
  
  # Create structured content for dynamic reports
  content <- list(
    section_title = "GO Enrichment Analysis",
    analysis_date = as.character(analysis_date),
    summary_text = paste(summary_text, collapse = "\n"),
    
    # Data for potential charts/tables
    ontology_stats = ontology_summary,
    enriched_terms = enriched_terms,
    total_enriched = nrow(enriched_terms),
    contributing_genes = gene_count,
    
    # Key metrics for quick access
    metrics = list(
      bp_enriched = sum(ontology_summary$ontology == "BP" & ontology_summary$trending > 0),
      mf_enriched = sum(ontology_summary$ontology == "MF" & ontology_summary$trending > 0),
      cc_enriched = sum(ontology_summary$ontology == "CC" & ontology_summary$trending > 0),
      highly_significant = if(nrow(enriched_terms) > 0) sum(enriched_terms$p_adjusted < 0.01) else 0,
      significant = if(nrow(enriched_terms) > 0) sum(enriched_terms$p_adjusted < 0.05) else 0,
      trending = if(nrow(enriched_terms) > 0) sum(enriched_terms$p_adjusted < 0.1) else 0
    )
  )
  
  return(content)
}