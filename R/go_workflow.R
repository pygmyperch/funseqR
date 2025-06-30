#' GO Enrichment Workflow Functions
#'
#' High-level wrapper functions for complete GO enrichment workflows
#'

#' Over-Representation Analysis (ORA) workflow for GO and KEGG enrichment
#'
#' @param con Database connection object
#' @param candidate_vcf_file Character. Path to candidate VCF file
#' @param background_file_id Integer. File ID of background dataset (or NULL to auto-detect)
#' @param blast_param_id Integer. Optional. Specific BLAST run ID to use for both datasets.
#'   If NULL, uses all available annotations. Default is NULL.
#' @param annotation_type Character. Type of annotations to analyze: "GO", "KEGG", "Pfam", "InterPro", "eggNOG", "both", or "all". Default is "both"
#' @param ontologies Character vector. GO ontologies to test: c("BP", "MF", "CC"). Only used when GO is included. Default is c("BP", "MF", "CC")
#' @param min_genes Integer. Minimum genes for term/pathway testing. Default is 5
#' @param max_genes Integer. Maximum genes for term/pathway testing. Default is 500
#' @param significance_threshold Numeric. FDR threshold for significance. Default is 0.05
#' @param method Character. Enrichment method: "clusterprofiler" or "legacy". Default is "clusterprofiler"
#' @param store_results Logical. Store results in database. Default is TRUE
#' @param create_plots Logical. Generate visualization plots. Default is TRUE
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return List containing all analysis results, plots, locus information, and summary tables:
#' \\itemize{
#'   \\item status: Analysis completion status
#'   \\item summary: Workflow summary with parameters and statistics
#'   \\item candidate_import: Information about imported candidate file
#'   \\item annotation_data: Raw annotation data used for enrichment
#'   \\item enrichment_results: Enrichment test results for each annotation type
#'   \\item locus_info: Detailed table of loci with annotations and dataset membership
#'   \\item plots: Generated visualization plots
#'   \\item analysis_ids: Database IDs of stored enrichment analyses
#' }
#'
#' @details
#' This function performs complete Over-Representation Analysis (ORA) workflow:
#' 1. Import candidate loci file
#' 2. Link to existing annotations
#' 3. Extract annotations for both datasets (optionally filtered by BLAST run)
#' 4. Perform enrichment testing using clusterProfiler:
#'    - GO: Tests specified ontologies (BP, MF, CC)
#'    - KEGG: Tests pathway enrichment
#'    - both: Runs both GO and KEGG analyses
#' 5. Create visualizations
#' 6. Store results in unified database schema
#' 
#' When blast_param_id is specified, only annotations from that specific BLAST run
#' are used for both candidate and background datasets, ensuring methodological 
#' consistency and enabling comparison of different annotation strategies.
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#' 
#' # Run both GO and KEGG enrichment (default)
#' results <- run_ORA(con, "candidates.vcf")
#' 
#' # Run only GO enrichment
#' go_results <- run_ORA(con, "candidates.vcf", annotation_type = "GO")
#' 
#' # Run only KEGG enrichment
#' kegg_results <- run_ORA(con, "candidates.vcf", annotation_type = "KEGG")
#' 
#' # Use custom significance threshold
#' results_lenient <- run_ORA(con, "candidates.vcf", significance_threshold = 0.1)
#' 
#' # Use only ORF-based annotations for both GO and KEGG
#' results_orf <- run_ORA(con, "candidates.vcf", blast_param_id = 1)
#' 
#' print(results$summary)
#' print(results$plots$GO_BP_bubble)
#' print(results$plots$KEGG_pathway_bubble)
#' }
#'
#' @export
run_ORA <- function(con, candidate_vcf_file, background_file_id = NULL,
                   blast_param_id = NULL, annotation_type = c("both", "GO", "KEGG", "Pfam", "InterPro", "eggNOG", "all"),
                   ontologies = c("BP", "MF", "CC"), 
                   min_genes = 5, max_genes = 500, significance_threshold = 0.05,
                   method = "clusterprofiler", store_results = TRUE, create_plots = TRUE, verbose = TRUE) {
  
  # Validate annotation_type parameter
  annotation_type <- match.arg(annotation_type)
  
  if (verbose) message("=== Starting Over-Representation Analysis (ORA) ===")
  if (verbose) message("Annotation types: ", annotation_type)
  
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
  
  # Step 2: Extract annotation data and perform enrichment
  if (verbose) message("\n=== Step 2: Extracting Annotation Data ===")
  
  enrichment_results <- list()
  enrichment_ids <- list()
  annotation_data <- list()
  
  # Determine which analyses to run
  run_go <- annotation_type %in% c("GO", "both", "all")
  run_kegg <- annotation_type %in% c("KEGG", "both", "all")
  run_pfam <- annotation_type %in% c("Pfam", "all")
  run_interpro <- annotation_type %in% c("InterPro", "all")
  run_eggnog <- annotation_type %in% c("eggNOG", "all")
  
  if (run_go) {
    if (verbose) message("  - Extracting GO terms...")
    go_data <- extract_go_terms_for_enrichment(con, candidate_import$file_id, background_file_id, 
                                             blast_param_id = blast_param_id, verbose = verbose)
    annotation_data[["GO"]] <- go_data
    
    if (length(go_data$foreground$genes) == 0) {
      warning("No GO terms found for candidate genes.")
      enrichment_results[["GO"]] <- list()
    } else {
      # Perform GO enrichment for each ontology
      if (verbose) message("\n=== Step 3a: Performing GO Enrichment Analysis ===")
      go_results <- list()
      go_ids <- list()
      
      for (ontology in ontologies) {
        if (verbose) message("  - Analyzing ", ontology, " ontology...")
        
        results <- perform_go_enrichment(go_data, ontology, min_genes = min_genes, 
                                       max_genes = max_genes, significance_threshold = significance_threshold,
                                       method = method, verbose = verbose)
        
        go_results[[ontology]] <- results
        
        # Store results in database if requested
        if (store_results && nrow(results) > 0) {
          analysis_id <- store_ora_results(
            con, candidate_import$file_id, background_file_id,
            results, "GO", ontology, 
            parameters = list(min_genes = min_genes, max_genes = max_genes, significance_threshold = significance_threshold),
            method = method,
            blast_param_id = blast_param_id,
            verbose = verbose
          )
          go_ids[[ontology]] <- analysis_id
        }
      }
      
      enrichment_results[["GO"]] <- go_results
      enrichment_ids[["GO"]] <- go_ids
    }
  }
  
  if (run_kegg) {
    if (verbose) message("  - Extracting KEGG pathways...")
    kegg_data <- extract_kegg_terms_for_enrichment(con, candidate_import$file_id, background_file_id, 
                                                 blast_param_id = blast_param_id, verbose = verbose)
    annotation_data[["KEGG"]] <- kegg_data
    
    if (length(kegg_data$foreground$genes) == 0) {
      warning("No KEGG pathways found for candidate genes.")
      enrichment_results[["KEGG"]] <- list()
    } else {
      # Perform KEGG enrichment
      if (verbose) message("\n=== Step 3b: Performing KEGG Enrichment Analysis ===")
      
      kegg_results <- perform_kegg_enrichment(kegg_data, min_genes = min_genes, 
                                            max_genes = max_genes, significance_threshold = significance_threshold,
                                            method = method, verbose = verbose)
      
      enrichment_results[["KEGG"]] <- list(PATHWAY = kegg_results)
      
      # Store results in database if requested
      if (store_results && nrow(kegg_results) > 0) {
        analysis_id <- store_ora_results(
          con, candidate_import$file_id, background_file_id,
          kegg_results, "KEGG", "PATHWAY", 
          parameters = list(min_genes = min_genes, max_genes = max_genes, significance_threshold = significance_threshold),
          method = method,
          blast_param_id = blast_param_id,
          verbose = verbose
        )
        enrichment_ids[["KEGG"]] <- list(PATHWAY = analysis_id)
      }
    }
  }
  
  # Step 4: Create visualizations
  plots <- list()
  if (create_plots) {
    if (verbose) message("\n=== Step 4: Creating Visualizations ===")
    
    # Create GO plots
    if (run_go && "GO" %in% names(enrichment_results)) {
      go_results <- enrichment_results[["GO"]]
      
      for (ontology in names(go_results)) {
        results <- go_results[[ontology]]
        
        if (nrow(results) > 0) {
          if (verbose) message("  - Creating plots for GO ", ontology, " ontology...")
          
          # Bubble plot
          plots[[paste0("GO_", ontology, "_bubble")]] <- create_go_bubble_plot(results)
          
          # Treemap (if treemapify is available)
          if (requireNamespace("treemapify", quietly = TRUE)) {
            plots[[paste0("GO_", ontology, "_treemap")]] <- create_go_treemap(results)
          }
          
          # Summary table
          plots[[paste0("GO_", ontology, "_table")]] <- create_go_summary_table(results)
        }
      }
      
      # Multi-ontology comparison plot if multiple ontologies tested
      if (length(go_results) > 1) {
        bp_res <- if ("BP" %in% names(go_results)) go_results[["BP"]] else NULL
        mf_res <- if ("MF" %in% names(go_results)) go_results[["MF"]] else NULL
        cc_res <- if ("CC" %in% names(go_results)) go_results[["CC"]] else NULL
        
        plots[["GO_comparison"]] <- create_go_comparison_plot(bp_res, mf_res, cc_res)
      }
    }
    
    # Create KEGG plots
    if (run_kegg && "KEGG" %in% names(enrichment_results)) {
      kegg_results <- enrichment_results[["KEGG"]][["PATHWAY"]]
      
      if (nrow(kegg_results) > 0) {
        if (verbose) message("  - Creating plots for KEGG pathways...")
        
        # Bubble plot (adapt GO plot for KEGG)
        plots[["KEGG_pathway_bubble"]] <- create_kegg_bubble_plot(kegg_results)
        
        # Summary table
        plots[["KEGG_pathway_table"]] <- create_kegg_summary_table(kegg_results)
      }
    }
  }
  
  # Step 5: Create summary
  if (verbose) message("\n=== Step 5: Creating Summary ===")
  
  # Build summary statistics for all analyses
  summary_stats_list <- list()
  
  if (run_go && "GO" %in% names(enrichment_results)) {
    go_results <- enrichment_results[["GO"]]
    for (ontology in names(go_results)) {
      results <- go_results[[ontology]]
      summary_stats_list[[paste0("GO_", ontology)]] <- data.frame(
        Analysis_Type = "GO",
        Term_Type = ontology,
        Terms_Tested = nrow(results),
        Significant_Terms = sum(results$p_adjusted < significance_threshold, na.rm = TRUE),
        Highly_Significant = sum(results$p_adjusted < (significance_threshold / 5), na.rm = TRUE),
        Top_Enrichment = if (nrow(results) > 0) round(max(results$fold_enrichment, na.rm = TRUE), 2) else 0,
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (run_kegg && "KEGG" %in% names(enrichment_results)) {
    kegg_results <- enrichment_results[["KEGG"]][["PATHWAY"]]
    summary_stats_list[["KEGG_PATHWAY"]] <- data.frame(
      Analysis_Type = "KEGG",
      Term_Type = "PATHWAY",
      Terms_Tested = nrow(kegg_results),
      Significant_Terms = sum(kegg_results$p_adjusted < significance_threshold, na.rm = TRUE),
      Highly_Significant = sum(kegg_results$p_adjusted < (significance_threshold / 5), na.rm = TRUE),
      Top_Enrichment = if (nrow(kegg_results) > 0) round(max(kegg_results$fold_enrichment, na.rm = TRUE), 2) else 0,
      stringsAsFactors = FALSE
    )
  }
  
  summary_stats <- do.call(rbind, summary_stats_list)
  
  # Calculate gene counts from annotation data
  total_foreground <- 0
  total_background <- 0
  
  if ("GO" %in% names(annotation_data)) {
    total_foreground <- max(total_foreground, length(annotation_data[["GO"]]$foreground$genes))
    total_background <- max(total_background, length(annotation_data[["GO"]]$background$genes))
  }
  
  if ("KEGG" %in% names(annotation_data)) {
    total_foreground <- max(total_foreground, length(annotation_data[["KEGG"]]$foreground$genes))
    total_background <- max(total_background, length(annotation_data[["KEGG"]]$background$genes))
  }
  
  workflow_summary <- list(
    analysis_date = Sys.time(),
    candidate_file = candidate_vcf_file,
    candidate_file_id = candidate_import$file_id,
    background_file_id = background_file_id,
    annotation_type = annotation_type,
    foreground_genes = total_foreground,
    background_genes = total_background,
    parameters = list(min_genes = min_genes, max_genes = max_genes, significance_threshold = significance_threshold),
    summary_stats = summary_stats,
    stored_analysis_ids = enrichment_ids
  )
  
  # Step 6: Create locus information table
  if (verbose) message("\n=== Step 6: Creating Locus Information Table ===")
  
  locus_info <- .create_locus_info_table(con, candidate_import$file_id, background_file_id, 
                                         blast_param_id, annotation_type, verbose)
  
  if (verbose) {
    message("=== ORA Workflow Complete ===")
    message("Annotation types analyzed: ", annotation_type)
    if (total_foreground > 0) {
      message("Candidate genes with annotations: ", total_foreground)
      message("Background genes with annotations: ", total_background)
    }
    if (!is.null(summary_stats)) {
      message("Total significant terms: ", sum(summary_stats$Significant_Terms, na.rm = TRUE))
    }
    if (store_results) {
      message("Results stored with analysis IDs: ", paste(unlist(enrichment_ids), collapse = ", "))
    }
    if (!is.null(locus_info)) {
      message("Locus information table created with ", nrow(locus_info), " entries")
    }
  }
  
  return(list(
    status = "success",
    summary = workflow_summary,
    candidate_import = candidate_import,
    annotation_data = annotation_data,
    enrichment_results = enrichment_results,
    locus_info = locus_info,
    plots = plots,
    analysis_ids = enrichment_ids
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

#' Filter locus information for significantly enriched terms
#'
#' Extracts loci that are associated with significantly enriched functional terms
#' from ORA (Over-Representation Analysis) results.
#'
#' @param ora_results List. Results from run_ORA() function
#' @param significance_threshold Numeric. FDR threshold for significance. Default is 0.05
#' @param annotation_types Character vector. Which annotation types to consider: "GO", "KEGG", "Pfam", "InterPro", "eggNOG", or "all". Default is "all"
#' @param dataset_type Character. Filter by dataset type: "candidate", "background", or "both". Default is "both"
#' @param return_summary Logical. Include summary statistics. Default is TRUE
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return List containing:
#' \itemize{
#'   \item enriched_loci: Data frame of loci associated with enriched terms
#'   \item enriched_terms: List of significantly enriched terms by annotation type
#'   \item summary: Summary statistics (if return_summary = TRUE)
#' }
#'
#' @details
#' This function processes ORA results to identify which genomic loci are
#' associated with significantly enriched functional terms. It works across
#' all annotation types and provides flexible filtering options.
#'
#' The function:
#' \enumerate{
#'   \item Extracts significantly enriched terms from each annotation type
#'   \item Maps these terms back to their associated genes/proteins
#'   \item Filters the locus_info table to include only enriched loci
#'   \item Provides summary statistics and enriched term lists
#' }
#'
#' @examples
#' \dontrun{
#' # Run ORA analysis
#' ora_results <- run_ORA(con, "candidates.vcf", annotation_type = "both")
#' 
#' # Get all enriched loci
#' enriched <- filter_enriched_loci(ora_results)
#' head(enriched$enriched_loci)
#' 
#' # Get only GO-enriched candidate loci
#' go_candidates <- filter_enriched_loci(
#'   ora_results, 
#'   annotation_types = "GO",
#'   dataset_type = "candidate",
#'   significance_threshold = 0.1
#' )
#' 
#' # View summary
#' print(enriched$summary)
#' }
#'
#' @export
filter_enriched_loci <- function(ora_results, significance_threshold = 0.05, 
                                 annotation_types = "all", dataset_type = "both",
                                 return_summary = TRUE, verbose = TRUE) {
  
  if (verbose) message("Filtering loci for significantly enriched terms...")
  
  # Validate inputs
  if (!"locus_info" %in% names(ora_results)) {
    stop("ORA results must contain locus_info table")
  }
  
  if (!"enrichment_results" %in% names(ora_results)) {
    stop("ORA results must contain enrichment_results")
  }
  
  locus_info <- ora_results$locus_info
  enrichment_results <- ora_results$enrichment_results
  
  # Determine which annotation types to process
  available_types <- names(enrichment_results)
  if (annotation_types[1] == "all") {
    annotation_types <- available_types
  } else {
    annotation_types <- intersect(annotation_types, available_types)
  }
  
  if (verbose) message("  - Processing annotation types: ", paste(annotation_types, collapse = ", "))
  
  # Extract enriched terms and associated genes
  enriched_genes <- c()
  enriched_terms_list <- list()
  
  for (ann_type in annotation_types) {
    if (verbose) message("    - Processing ", ann_type, " enrichment results...")
    
    type_results <- enrichment_results[[ann_type]]
    
    if (ann_type == "GO") {
      # GO has multiple ontologies
      for (ontology in names(type_results)) {
        results_df <- type_results[[ontology]]
        if (nrow(results_df) > 0) {
          significant <- results_df[results_df$p_adjusted < significance_threshold, ]
          if (nrow(significant) > 0) {
            # Extract gene IDs from gene_ids column
            genes <- unique(unlist(strsplit(significant$gene_ids, "/")))
            enriched_genes <- c(enriched_genes, genes)
            enriched_terms_list[[paste0(ann_type, "_", ontology)]] <- significant
            if (verbose) message("      - Found ", nrow(significant), " enriched ", ontology, " terms")
          }
        }
      }
    } else {
      # Other annotation types (KEGG, Pfam, etc.)
      if (is.list(type_results) && length(type_results) > 0) {
        # Handle case where results are nested (like KEGG$PATHWAY)
        results_df <- if ("PATHWAY" %in% names(type_results)) {
          type_results$PATHWAY
        } else {
          type_results[[1]]
        }
        
        if (is.data.frame(results_df) && nrow(results_df) > 0) {
          significant <- results_df[results_df$p_adjusted < significance_threshold, ]
          if (nrow(significant) > 0) {
            genes <- unique(unlist(strsplit(significant$gene_ids, "/")))
            enriched_genes <- c(enriched_genes, genes)
            enriched_terms_list[[ann_type]] <- significant
            if (verbose) message("      - Found ", nrow(significant), " enriched ", ann_type, " terms")
          }
        }
      }
    }
  }
  
  # Remove duplicates and clean gene IDs
  enriched_genes <- unique(enriched_genes)
  enriched_genes <- enriched_genes[!is.na(enriched_genes) & enriched_genes != ""]
  
  if (verbose) message("  - Total unique genes in enriched terms: ", length(enriched_genes))
  
  if (length(enriched_genes) == 0) {
    warning("No significantly enriched terms found")
    return(list(
      enriched_loci = data.frame(),
      enriched_terms = enriched_terms_list,
      summary = data.frame()
    ))
  }
  
  # Filter locus_info for enriched genes
  enriched_loci <- .filter_loci_by_genes(locus_info, enriched_genes, verbose)
  
  # Filter by dataset type if requested
  if (dataset_type != "both") {
    enriched_loci <- enriched_loci[enriched_loci$dataset_type == dataset_type, ]
    if (verbose) message("  - Filtered to ", dataset_type, " loci: ", nrow(enriched_loci))
  }
  
  # Create summary statistics
  summary_stats <- NULL
  if (return_summary) {
    summary_stats <- .create_enriched_loci_summary(enriched_loci, enriched_terms_list, locus_info)
  }
  
  if (verbose) {
    message("Filtering complete:")
    message("  - Enriched loci found: ", nrow(enriched_loci))
    message("  - Candidate loci: ", sum(enriched_loci$dataset_type == "candidate", na.rm = TRUE))
    message("  - Background loci: ", sum(enriched_loci$dataset_type == "background", na.rm = TRUE))
  }
  
  return(list(
    enriched_loci = enriched_loci,
    enriched_terms = enriched_terms_list,
    summary = summary_stats
  ))
}

#' Filter loci by gene IDs
#' @keywords internal
.filter_loci_by_genes <- function(locus_info, enriched_genes, verbose) {
  
  if (nrow(locus_info) == 0 || length(enriched_genes) == 0) {
    return(data.frame())
  }
  
  # Create a logical vector for matching loci
  matches <- rep(FALSE, nrow(locus_info))
  
  for (i in 1:nrow(locus_info)) {
    # Check uniprot_accessions column
    if (!is.na(locus_info$uniprot_accessions[i])) {
      accessions <- strsplit(locus_info$uniprot_accessions[i], ";")[[1]]
      if (any(accessions %in% enriched_genes)) {
        matches[i] <- TRUE
        next
      }
    }
    
    # Check gene_names column if available
    if ("gene_names" %in% names(locus_info) && !is.na(locus_info$gene_names[i])) {
      gene_names <- strsplit(locus_info$gene_names[i], ";")[[1]]
      if (any(gene_names %in% enriched_genes)) {
        matches[i] <- TRUE
      }
    }
  }
  
  filtered_loci <- locus_info[matches, ]
  
  if (verbose) message("  - Matched ", nrow(filtered_loci), " loci to enriched genes")
  
  return(filtered_loci)
}

#' Create summary statistics for enriched loci
#' @keywords internal
.create_enriched_loci_summary <- function(enriched_loci, enriched_terms_list, all_loci) {
  
  if (nrow(enriched_loci) == 0) {
    return(data.frame())
  }
  
  # Basic statistics
  total_loci <- nrow(all_loci)
  enriched_count <- nrow(enriched_loci)
  enriched_percentage <- round((enriched_count / total_loci) * 100, 2)
  
  candidate_enriched <- sum(enriched_loci$dataset_type == "candidate", na.rm = TRUE)
  background_enriched <- sum(enriched_loci$dataset_type == "background", na.rm = TRUE)
  
  total_candidates <- sum(all_loci$dataset_type == "candidate", na.rm = TRUE)
  total_background <- sum(all_loci$dataset_type == "background", na.rm = TRUE)
  
  candidate_enriched_pct <- if (total_candidates > 0) {
    round((candidate_enriched / total_candidates) * 100, 2)
  } else 0
  
  background_enriched_pct <- if (total_background > 0) {
    round((background_enriched / total_background) * 100, 2)
  } else 0
  
  # Count enriched terms by type
  term_counts <- data.frame(
    annotation_type = names(enriched_terms_list),
    enriched_terms = sapply(enriched_terms_list, nrow),
    stringsAsFactors = FALSE
  )
  
  # Create overall summary
  overall_summary <- data.frame(
    metric = c("Total loci", "Enriched loci", "Enriched percentage", 
               "Candidate loci (total)", "Candidate loci (enriched)", "Candidate enriched %",
               "Background loci (total)", "Background loci (enriched)", "Background enriched %"),
    value = c(total_loci, enriched_count, enriched_percentage,
              total_candidates, candidate_enriched, candidate_enriched_pct,
              total_background, background_enriched, background_enriched_pct),
    stringsAsFactors = FALSE
  )
  
  return(list(
    overall = overall_summary,
    by_annotation_type = term_counts
  ))
}

#' Create locus information table for ORA results
#' @keywords internal
.create_locus_info_table <- function(con, candidate_file_id, background_file_id, blast_param_id, annotation_type, verbose) {
  
  if (verbose) message("  - Creating comprehensive locus information table...")
  
  # Build base query conditions
  where_conditions <- c()
  params <- list()
  
  if (!is.null(blast_param_id)) {
    where_conditions <- c(where_conditions, "br.blast_param_id = ?")
    params <- c(params, list(blast_param_id))
  }
  
  where_clause <- if (length(where_conditions) > 0) {
    paste("AND", paste(where_conditions, collapse = " AND "))
  } else {
    ""
  }
  
  # Get candidate and background loci information with annotations
  locus_query <- paste0(
    "SELECT DISTINCT ",
    "vd.vcf_id || '_' || vd.chromosome || '_' || vd.position as locus_id, ",
    "vd.chromosome, ",
    "vd.position, ",
    "vd.ref, ",
    "vd.alt, ",
    "CASE ",
    "  WHEN vd.file_id = ? THEN 'candidate' ",
    "  WHEN vd.file_id = ? THEN 'background' ",
    "  ELSE 'other' ",
    "END as dataset_type, ",
    "a.uniprot_accession, ",
    "a.entry_name, ",
    "a.gene_names, ",
    "br.e_value, ",
    "br.bit_score, ",
    "br.percent_identity ",
    "FROM vcf_data vd ",
    "JOIN flanking_sequences fs ON vd.vcf_id = fs.vcf_id ",
    "JOIN blast_results br ON fs.flanking_id = br.flanking_id ",
    "JOIN annotations a ON br.blast_result_id = a.blast_result_id ",
    "WHERE (vd.file_id = ? OR vd.file_id = ?) ",
    where_clause, " ",
    "ORDER BY vd.chromosome, vd.position"
  )
  
  # Parameters: candidate_file_id, background_file_id, candidate_file_id, background_file_id, [blast_param_id]
  query_params <- c(list(candidate_file_id, background_file_id, candidate_file_id, background_file_id), params)
  
  # Execute query with error handling
  locus_data <- tryCatch({
    DBI::dbGetQuery(con, locus_query, query_params)
  }, error = function(e) {
    if (verbose) message("    - Error executing locus query: ", e$message)
    return(data.frame())
  })
  
  if (nrow(locus_data) == 0) {
    if (verbose) message("    - No loci found with annotations")
    return(data.frame())
  }
  
  if (verbose) message("    - Found ", nrow(locus_data), " annotated loci")
  
  # Check which annotation tables exist
  available_tables <- DBI::dbListTables(con)
  
  # Add annotation details based on annotation_type and table availability
  if (annotation_type %in% c("GO", "both", "all") && "go_terms" %in% available_tables) {
    locus_data <- .add_go_annotations_to_locus_table(con, locus_data, where_clause, params, verbose)
  }
  
  if (annotation_type %in% c("KEGG", "both", "all") && "kegg_references" %in% available_tables) {
    locus_data <- .add_kegg_annotations_to_locus_table(con, locus_data, where_clause, params, verbose)
  }
  
  if (annotation_type %in% c("Pfam", "all") && "pfam_domains" %in% available_tables) {
    locus_data <- .add_pfam_annotations_to_locus_table(con, locus_data, where_clause, params, verbose)
  } else if (annotation_type %in% c("Pfam", "all") && verbose) {
    message("    - Pfam table not found, skipping Pfam annotations")
  }
  
  if (annotation_type %in% c("InterPro", "all") && "interpro_families" %in% available_tables) {
    locus_data <- .add_interpro_annotations_to_locus_table(con, locus_data, where_clause, params, verbose)
  } else if (annotation_type %in% c("InterPro", "all") && verbose) {
    message("    - InterPro table not found, skipping InterPro annotations")
  }
  
  if (annotation_type %in% c("eggNOG", "all") && "eggnog_categories" %in% available_tables) {
    locus_data <- .add_eggnog_annotations_to_locus_table(con, locus_data, where_clause, params, verbose)
  } else if (annotation_type %in% c("eggNOG", "all") && verbose) {
    message("    - eggNOG table not found, skipping eggNOG annotations")
  }
  
  # Create summary by locus (aggregate multiple annotations per locus)
  summary_data <- tryCatch({
    locus_data %>%
      dplyr::group_by(locus_id, chromosome, position, ref, alt, dataset_type) %>%
      dplyr::summarise(
        uniprot_accessions = paste(unique(uniprot_accession), collapse = ";"),
        gene_names = paste(unique(na.omit(gene_names)), collapse = ";"),
        entry_names = paste(unique(na.omit(entry_name)), collapse = ";"),
        best_e_value = min(e_value, na.rm = TRUE),
        best_bit_score = max(bit_score, na.rm = TRUE),
        avg_percent_identity = round(mean(percent_identity, na.rm = TRUE), 2),
        annotation_count = dplyr::n(),
        .groups = "drop"
      )
  }, error = function(e) {
    if (verbose) message("    - Error processing locus data: ", e$message)
    return(data.frame())
  })
  
  # Add aggregated annotation information if present
  if ("go_terms" %in% names(locus_data)) {
    go_summary <- locus_data %>%
      dplyr::group_by(locus_id) %>%
      dplyr::summarise(
        go_terms = paste(unique(na.omit(unlist(strsplit(go_terms, ";")))), collapse = ";"),
        go_categories = paste(unique(na.omit(unlist(strsplit(go_categories, ";")))), collapse = ";"),
        .groups = "drop"
      )
    
    summary_data <- merge(summary_data, go_summary, by = "locus_id", all.x = TRUE)
  }
  
  if ("kegg_pathways" %in% names(locus_data)) {
    kegg_summary <- locus_data %>%
      dplyr::group_by(locus_id) %>%
      dplyr::summarise(
        kegg_pathways = paste(unique(na.omit(unlist(strsplit(kegg_pathways, ";")))), collapse = ";"),
        .groups = "drop"
      )
    
    summary_data <- merge(summary_data, kegg_summary, by = "locus_id", all.x = TRUE)
  }
  
  if (verbose) message("    - Created locus table with ", nrow(summary_data), " unique loci")
  
  return(as.data.frame(summary_data))
}

#' Add GO annotations to locus table
#' @keywords internal
.add_go_annotations_to_locus_table <- function(con, locus_data, where_clause, params, verbose) {
  
  # Get GO terms for the loci
  locus_ids_str <- paste0("'", paste(unique(locus_data$locus_id), collapse = "', '"), "'")
  
  go_query <- paste0(
    "SELECT DISTINCT ",
    "vd.vcf_id || '_' || vd.chromosome || '_' || vd.position as locus_id, ",
    "gt.go_id, ",
    "gt.go_term, ",
    "gt.go_category ",
    "FROM vcf_data vd ",
    "JOIN flanking_sequences fs ON vd.vcf_id = fs.vcf_id ",
    "JOIN blast_results br ON fs.flanking_id = br.flanking_id ",
    "JOIN annotations a ON br.blast_result_id = a.blast_result_id ",
    "JOIN go_terms gt ON a.annotation_id = gt.annotation_id ",
    "WHERE vd.vcf_id || '_' || vd.chromosome || '_' || vd.position IN (", locus_ids_str, ") ",
    if (nchar(where_clause) > 0) paste(" ", where_clause) else ""
  )
  
  go_data <- if (length(params) > 0) {
    DBI::dbGetQuery(con, go_query, params)
  } else {
    DBI::dbGetQuery(con, go_query)
  }
  
  if (nrow(go_data) > 0) {
    # Aggregate GO terms by locus
    go_summary <- go_data %>%
      dplyr::group_by(locus_id) %>%
      dplyr::summarise(
        go_terms = paste(unique(go_id), collapse = ";"),
        go_names = paste(unique(go_term), collapse = ";"),
        go_categories = paste(unique(go_category), collapse = ";"),
        .groups = "drop"
      )
    
    # Merge with locus data
    locus_data <- merge(locus_data, go_summary, by = "locus_id", all.x = TRUE)
  } else {
    locus_data$go_terms <- NA_character_
    locus_data$go_names <- NA_character_
    locus_data$go_categories <- NA_character_
  }
  
  return(locus_data)
}

#' Add KEGG annotations to locus table
#' @keywords internal
.add_kegg_annotations_to_locus_table <- function(con, locus_data, where_clause, params, verbose) {
  
  # Get KEGG pathways for the loci
  locus_ids_str <- paste0("'", paste(unique(locus_data$locus_id), collapse = "', '"), "'")
  
  kegg_query <- paste0(
    "SELECT DISTINCT ",
    "vd.vcf_id || '_' || vd.chromosome || '_' || vd.position as locus_id, ",
    "kr.kegg_id, ",
    "kr.pathway_name ",
    "FROM vcf_data vd ",
    "JOIN flanking_sequences fs ON vd.vcf_id = fs.vcf_id ",
    "JOIN blast_results br ON fs.flanking_id = br.flanking_id ",
    "JOIN annotations a ON br.blast_result_id = a.blast_result_id ",
    "JOIN kegg_references kr ON a.annotation_id = kr.annotation_id ",
    "WHERE vd.vcf_id || '_' || vd.chromosome || '_' || vd.position IN (", locus_ids_str, ") ",
    if (nchar(where_clause) > 0) paste(" ", where_clause) else ""
  )
  
  kegg_data <- if (length(params) > 0) {
    DBI::dbGetQuery(con, kegg_query, params)
  } else {
    DBI::dbGetQuery(con, kegg_query)
  }
  
  if (nrow(kegg_data) > 0) {
    # Aggregate KEGG pathways by locus
    kegg_summary <- kegg_data %>%
      dplyr::group_by(locus_id) %>%
      dplyr::summarise(
        kegg_pathways = paste(unique(kegg_id), collapse = ";"),
        kegg_pathway_names = paste(unique(pathway_name), collapse = ";"),
        .groups = "drop"
      )
    
    # Merge with locus data
    locus_data <- merge(locus_data, kegg_summary, by = "locus_id", all.x = TRUE)
  } else {
    locus_data$kegg_pathways <- NA_character_
    locus_data$kegg_pathway_names <- NA_character_
  }
  
  return(locus_data)
}

#' Add Pfam annotations to locus table (placeholder)
#' @keywords internal
.add_pfam_annotations_to_locus_table <- function(con, locus_data, where_clause, params, verbose) {
  # Placeholder for Pfam annotations - similar structure to GO/KEGG
  locus_data$pfam_domains <- NA_character_
  return(locus_data)
}

#' Add InterPro annotations to locus table (placeholder)
#' @keywords internal
.add_interpro_annotations_to_locus_table <- function(con, locus_data, where_clause, params, verbose) {
  # Placeholder for InterPro annotations - similar structure to GO/KEGG
  locus_data$interpro_families <- NA_character_
  return(locus_data)
}

#' Add eggNOG annotations to locus table (placeholder)
#' @keywords internal
.add_eggnog_annotations_to_locus_table <- function(con, locus_data, where_clause, params, verbose) {
  # Placeholder for eggNOG annotations - similar structure to GO/KEGG
  locus_data$eggnog_categories <- NA_character_
  return(locus_data)
}