#' GO Enrichment Workflow Functions
#'
#' High-level wrapper functions for complete GO enrichment workflows
#'

#' Complete GO enrichment analysis workflow
#'
#' @param con Database connection object
#' @param candidate_vcf_file Character. Path to candidate VCF file
#' @param background_file_id Integer. File ID of background dataset (or NULL to auto-detect)
#' @param blast_param_id Integer. Optional. Specific BLAST run ID to use for both datasets.
#'   If NULL, uses all available annotations. Default is NULL.
#' @param ontologies Character vector. GO ontologies to test: c("BP", "MF", "CC"). Default is c("BP", "MF", "CC")
#' @param min_genes Integer. Minimum genes for GO term testing. Default is 5
#' @param max_genes Integer. Maximum genes for GO term testing. Default is 500
#' @param significance_threshold Numeric. FDR threshold for significance. Default is 0.05
#' @param store_results Logical. Store results in database. Default is TRUE
#' @param create_plots Logical. Generate visualization plots. Default is TRUE
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return List containing all analysis results, plots, and summary tables
#'
#' @details
#' This function performs a complete GO enrichment workflow:
#' 1. Import candidate loci file
#' 2. Link to existing annotations
#' 3. Extract GO terms for both datasets (optionally filtered by BLAST run)
#' 4. Perform enrichment testing for specified ontologies
#' 5. Create visualizations
#' 6. Store results in database (optional)
#' 
#' When blast_param_id is specified, only annotations from that specific BLAST run
#' are used for both candidate and background datasets, ensuring methodological 
#' consistency and enabling comparison of different annotation strategies.
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#' 
#' # Use all available annotations and all ontologies (default behavior)
#' results <- run_go_enrichment_workflow(con, "candidates.vcf")
#' 
#' # Use custom significance threshold (e.g., 0.1 for more lenient)
#' results_lenient <- run_go_enrichment_workflow(con, "candidates.vcf", significance_threshold = 0.1)
#' 
#' # Use only specific ontologies and stricter threshold
#' results_strict <- run_go_enrichment_workflow(con, "candidates.vcf", 
#'                                             ontologies = c("BP", "MF"), 
#'                                             significance_threshold = 0.01)
#' 
#' # Use only ORF-based annotations
#' results_orf <- run_go_enrichment_workflow(con, "candidates.vcf", blast_param_id = 1)
#' 
#' print(results_lenient$summary)
#' print(results_lenient$plots$BP_bubble)
#' }
#'
#' @export
run_go_enrichment_workflow <- function(con, candidate_vcf_file, background_file_id = NULL,
                                      blast_param_id = NULL, ontologies = c("BP", "MF", "CC"), 
                                      min_genes = 5, max_genes = 500, significance_threshold = 0.05,
                                      method = "clusterprofiler", store_results = TRUE, create_plots = TRUE, verbose = TRUE) {
  
  if (verbose) message("=== Starting GO Enrichment Workflow ===")
  
  # Report BLAST run configuration
  if (!is.null(blast_param_id)) {
    if (verbose) {
      message("BLAST run filtering:")
      message("  - Using BLAST run ID: ", blast_param_id)
    }
  } else {
    if (verbose) {
      message("Using all available annotations")
    }
  }
  
  # Auto-detect background file if not provided
  if (is.null(background_file_id)) {
    if (verbose) message("Auto-detecting background dataset...")
    
    # Get all VCF files excluding the candidate file
    input_files <- DBI::dbGetQuery(con, 
      "SELECT file_id, file_name FROM input_files WHERE file_type = 'vcf'")
    
    # Exclude files with "candidate" or similar in the name
    candidate_patterns <- c("candidate", "adaptive", "outlier", "fst", "selection")
    background_files <- input_files[!grepl(paste(candidate_patterns, collapse = "|"), 
                                          tolower(input_files$file_name)), ]
    
    if (nrow(background_files) == 0) {
      stop("No suitable background dataset found. Please specify background_file_id.")
    }
    
    background_file_id <- background_files$file_id[1]
    if (verbose) message("  - Using background file: ", background_files$file_name[1], " (ID: ", background_file_id, ")")
  }
  
  # Step 1: Import candidate loci
  if (verbose) message("\n=== Step 1: Importing Candidate Loci ===")
  candidate_import <- import_candidate_loci(con, candidate_vcf_file, background_file_id, verbose = verbose)
  
  # Step 2: Extract GO terms
  if (verbose) message("\n=== Step 2: Extracting GO Terms ===")
  go_data <- extract_go_terms_for_enrichment(con, candidate_import$file_id, background_file_id, 
                                           blast_param_id = blast_param_id, verbose = verbose)
  
  if (length(go_data$foreground$genes) == 0) {
    warning("No GO terms found for candidate genes. Cannot perform enrichment analysis.")
    return(list(
      status = "failed",
      message = "No GO terms found for candidate genes",
      candidate_import = candidate_import
    ))
  }
  
  # Step 3: Perform enrichment analysis for each ontology
  if (verbose) message("\n=== Step 3: Performing Enrichment Analysis ===")
  enrichment_results <- list()
  enrichment_ids <- list()
  
  for (ontology in ontologies) {
    if (verbose) message("  - Analyzing ", ontology, " ontology...")
    
    results <- perform_go_enrichment(go_data, ontology, min_genes = min_genes, 
                                   max_genes = max_genes, significance_threshold = significance_threshold,
                                   method = method, verbose = verbose)
    
    enrichment_results[[ontology]] <- results
    
    # Store results in database if requested
    if (store_results && nrow(results) > 0) {
      enrichment_id <- store_go_enrichment_results(
        con, candidate_import$file_id, background_file_id,
        results, ontology, 
        parameters = list(min_genes = min_genes, max_genes = max_genes, significance_threshold = significance_threshold),
        method = method,
        blast_param_id = blast_param_id,
        verbose = verbose
      )
      enrichment_ids[[ontology]] <- enrichment_id
    }
  }
  
  # Step 4: Create visualizations
  plots <- list()
  if (create_plots) {
    if (verbose) message("\n=== Step 4: Creating Visualizations ===")
    
    for (ontology in names(enrichment_results)) {
      results <- enrichment_results[[ontology]]
      
      if (nrow(results) > 0) {
        if (verbose) message("  - Creating plots for ", ontology, " ontology...")
        
        # Bubble plot
        plots[[paste0(ontology, "_bubble")]] <- create_go_bubble_plot(results)
        
        # Treemap (if treemapify is available)
        if (requireNamespace("treemapify", quietly = TRUE)) {
          plots[[paste0(ontology, "_treemap")]] <- create_go_treemap(results)
        }
        
        # Summary table
        plots[[paste0(ontology, "_table")]] <- create_go_summary_table(results)
      }
    }
    
    # Multi-ontology comparison plot if multiple ontologies tested
    if (length(enrichment_results) > 1) {
      bp_res <- if ("BP" %in% names(enrichment_results)) enrichment_results[["BP"]] else NULL
      mf_res <- if ("MF" %in% names(enrichment_results)) enrichment_results[["MF"]] else NULL
      cc_res <- if ("CC" %in% names(enrichment_results)) enrichment_results[["CC"]] else NULL
      
      plots[["comparison"]] <- create_go_comparison_plot(bp_res, mf_res, cc_res)
    }
  }
  
  # Step 5: Create summary
  if (verbose) message("\n=== Step 5: Creating Summary ===")
  
  summary_stats <- data.frame(
    Ontology = names(enrichment_results),
    Terms_Tested = sapply(enrichment_results, nrow),
    Significant_Terms = sapply(enrichment_results, function(x) sum(x$p_adjusted < significance_threshold, na.rm = TRUE)),
    Highly_Significant = sapply(enrichment_results, function(x) sum(x$p_adjusted < (significance_threshold / 5), na.rm = TRUE)),
    Top_Enrichment = sapply(enrichment_results, function(x) {
      if (nrow(x) > 0) round(max(x$fold_enrichment, na.rm = TRUE), 2) else 0
    }),
    stringsAsFactors = FALSE
  )
  
  workflow_summary <- list(
    analysis_date = Sys.time(),
    candidate_file = candidate_vcf_file,
    candidate_file_id = candidate_import$file_id,
    background_file_id = background_file_id,
    foreground_genes = length(go_data$foreground$genes),
    background_genes = length(go_data$background$genes),
    ontologies_tested = ontologies,
    parameters = list(min_genes = min_genes, max_genes = max_genes, significance_threshold = significance_threshold),
    summary_stats = summary_stats,
    stored_analysis_ids = enrichment_ids
  )
  
  if (verbose) {
    message("=== Workflow Complete ===")
    message("Candidate genes with annotations: ", workflow_summary$foreground_genes)
    message("Background genes with annotations: ", workflow_summary$background_genes)
    message("Total significant terms: ", sum(summary_stats$Significant_Terms))
    if (store_results) {
      message("Results stored with IDs: ", paste(unlist(enrichment_ids), collapse = ", "))
    }
  }
  
  return(list(
    status = "success",
    summary = workflow_summary,
    candidate_import = candidate_import,
    go_data = go_data,
    enrichment_results = enrichment_results,
    plots = plots,
    enrichment_ids = enrichment_ids
  ))
}

#' Generate GO enrichment report section
#'
#' @param con Database connection object
#' @param project_id Integer. Project ID
#' @param candidate_file_pattern Character. Pattern to identify candidate files. Default is "candidate|adaptive|outlier"
#' @param max_terms_plot Integer. Maximum terms to show in plots. Default is 15
#' @param include_treemap Logical. Include treemap plots. Default is TRUE
#'
#' @return Character vector of R Markdown code for inclusion in reports
#'
#' @details
#' Generates R Markdown code that can be included in analysis reports.
#' Automatically detects candidate and background files and runs the analysis.
#'
#' @export
generate_go_enrichment_report_section <- function(con, project_id, 
                                                  candidate_file_pattern = "candidate|adaptive|outlier",
                                                  max_terms_plot = 15, include_treemap = TRUE) {
  
  # Create R Markdown code for GO enrichment section
  rmd_code <- c(
    "## GO Enrichment Analysis",
    "",
    "```{r go-enrichment-setup, include=FALSE}",
    "# Detect candidate and background files",
    "input_files <- DBI::dbGetQuery(con, ",
    "  'SELECT file_id, file_name FROM input_files WHERE project_id = ? AND file_type = \"vcf\"',",
    "  params = list(project_id))",
    "",
    "# Identify candidate files",
    paste0("candidate_files <- input_files[grepl('", candidate_file_pattern, "', tolower(input_files$file_name)), ]"),
    paste0("background_files <- input_files[!grepl('", candidate_file_pattern, "', tolower(input_files$file_name)), ]"),
    "",
    "go_enrichment_available <- nrow(candidate_files) > 0 && nrow(background_files) > 0",
    "```",
    "",
    "```{r go-enrichment-analysis, eval=go_enrichment_available}",
    "if (go_enrichment_available) {",
    "  # Extract GO terms for enrichment analysis",
    "  go_data <- extract_go_terms_for_enrichment(con, candidate_files$file_id[1], background_files$file_id[1], verbose = FALSE)",
    "  ",
    "  # Perform enrichment for each ontology",
    "  bp_results <- perform_go_enrichment(go_data, 'BP', verbose = FALSE)",
    "  mf_results <- perform_go_enrichment(go_data, 'MF', verbose = FALSE)",
    "  ",
    "  # Create summary statistics",
    "  enrichment_summary <- data.frame(",
    "    Ontology = c('Biological Process', 'Molecular Function'),",
    "    Candidate_Genes = rep(length(go_data$foreground$genes), 2),",
    "    Background_Genes = rep(length(go_data$background$genes), 2),",
    "    Terms_Tested = c(nrow(bp_results), nrow(mf_results)),",
    "    Significant = c(sum(bp_results$p_adjusted < 0.05, na.rm = TRUE), ",
    "                    sum(mf_results$p_adjusted < 0.05, na.rm = TRUE))",
    "  )",
    "  ",
    "  kable(enrichment_summary, caption = 'GO Enrichment Summary')",
    "}",
    "```",
    "",
    "### Biological Process Enrichment",
    "",
    "```{r bp-enrichment-plot, eval=go_enrichment_available, fig.height=8, fig.width=10}",
    "if (go_enrichment_available && nrow(bp_results) > 0) {",
    paste0("  bp_plot <- create_go_bubble_plot(bp_results, max_terms = ", max_terms_plot, ")"),
    "  print(bp_plot)",
    "} else {",
    "  cat('No significant Biological Process terms found.')",
    "}",
    "```",
    "",
    "### Molecular Function Enrichment", 
    "",
    "```{r mf-enrichment-plot, eval=go_enrichment_available, fig.height=6, fig.width=10}",
    "if (go_enrichment_available && nrow(mf_results) > 0) {",
    paste0("  mf_plot <- create_go_bubble_plot(mf_results, max_terms = ", max_terms_plot, ")"),
    "  print(mf_plot)",
    "} else {",
    "  cat('No significant Molecular Function terms found.')",
    "}",
    "```"
  )
  
  # Add treemap section if requested
  if (include_treemap) {
    treemap_code <- c(
      "",
      "### GO Term Treemaps",
      "",
      "```{r go-treemaps, eval=go_enrichment_available, fig.height=8, fig.width=12}",
      "if (go_enrichment_available && requireNamespace('treemapify', quietly = TRUE)) {",
      "  if (nrow(bp_results) > 0) {",
      "    bp_treemap <- create_go_treemap(bp_results)",
      "    print(bp_treemap)",
      "  }",
      "  ",
      "  if (nrow(mf_results) > 0) {",
      "    mf_treemap <- create_go_treemap(mf_results)",
      "    print(mf_treemap)",
      "  }",
      "} else {",
      "  cat('Treemap plots require the treemapify package.')",
      "}",
      "```"
    )
    rmd_code <- c(rmd_code, treemap_code)
  }
  
  # Add summary table
  summary_code <- c(
    "",
    "### Enrichment Results Table",
    "",
    "```{r go-summary-table, eval=go_enrichment_available}",
    "if (go_enrichment_available) {",
    "  # Combine top results from all ontologies",
    "  all_results <- rbind(",
    "    if (nrow(bp_results) > 0) bp_results[1:min(10, nrow(bp_results)), ] else NULL,",
    "    if (nrow(mf_results) > 0) mf_results[1:min(10, nrow(mf_results)), ] else NULL",
    "  )",
    "  ",
    "  if (!is.null(all_results) && nrow(all_results) > 0) {",
    "    summary_table <- create_go_summary_table(all_results, max_terms = 20)",
    "    kable(summary_table, caption = 'Top GO Enrichment Results')",
    "  } else {",
    "    cat('No significant enrichment results to display.')",
    "  }",
    "} else {",
    "  cat('GO enrichment analysis requires both candidate and background datasets.')",
    "}",
    "```"
  )
  
  rmd_code <- c(rmd_code, summary_code)
  
  return(rmd_code)
}