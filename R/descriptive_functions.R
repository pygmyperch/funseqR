#' Descriptive Functions for Functional Annotation Analysis
#'
#' This file contains functions to summarize and describe functional annotations
#' from BLAST results, GO terms, UniProt data, and enrichment analyses.
#' Enhanced with KEGGREST and ggkegg integration for improved pathway analysis.

#' Summarize GO annotations in the database
#'
#' Provides comprehensive statistics about GO term annotations across different
#' ontology categories (BP, MF, CC) and blast parameter sets, with optional filtering
#' for candidate loci and functional module grouping.
#'
#' @param con Database connection object
#' @param blast_param_id Integer. Specific BLAST parameter set to analyze. If NULL, analyzes all.
#' @param candidate_loci Data frame with chromosome and position columns, or NULL for all loci.
#' @param include_evidence Logical. Include GO evidence code statistics. Default is TRUE.
#' @param include_modules Logical. Group GO terms into functional modules. Default is FALSE.
#' @param min_frequency Integer. Minimum frequency threshold for GO terms to include. Default is 1.
#' @param verbose Logical. Print progress information. Default is TRUE.
#'
#' @return List containing GO annotation summaries:
#' \itemize{
#'   \item overview: Summary statistics by ontology category
#'   \item top_terms: Most frequent GO terms by category
#'   \item evidence_summary: GO evidence code distribution (if include_evidence = TRUE)
#'   \item coverage_stats: Coverage statistics across BLAST parameters
#'   \item functional_modules: GO terms grouped by functional modules (if include_modules = TRUE)
#' }
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#' 
#' # Analyze all GO annotations
#' go_summary <- summarize_go_annotations(con)
#' print(go_summary$overview)
#' 
#' # Analyze specific candidate loci with module grouping
#' candidates <- data.frame(chromosome = c("LG1", "LG2"), position = c(12345, 67890))
#' go_summary_candidates <- summarize_go_annotations(con, candidate_loci = candidates, 
#'                                                   include_modules = TRUE)
#' 
#' close_funseq_db(con)
#' }
#'
#' @export
summarize_go_annotations <- function(con, blast_param_id = NULL, candidate_loci = NULL,
                                   include_evidence = TRUE, include_modules = FALSE,
                                   min_frequency = 1, verbose = TRUE) {
  
  if (verbose) message("Summarizing GO annotations...")
  
  # Check if required tables exist
  tables <- DBI::dbListTables(con)
  required_tables <- c("go_terms", "annotations", "blast_results", "blast_parameters")
  missing_tables <- required_tables[!required_tables %in% tables]
  
  if (length(missing_tables) > 0) {
    stop("Missing required tables: ", paste(missing_tables, collapse = ", "))
  }
  
  # Build base query with optional filters
  base_conditions <- c()
  params <- list()
  
  if (!is.null(blast_param_id)) {
    base_conditions <- c(base_conditions, "bp.blast_param_id = ?")
    params <- c(params, list(blast_param_id))
    if (verbose) message("  - Filtering for blast_param_id: ", blast_param_id)
  }
  
  # Handle candidate loci filtering
  temp_table_name <- NULL
  loci_join <- ""
  
  if (!is.null(candidate_loci)) {
    if (!all(c("chromosome", "position") %in% colnames(candidate_loci))) {
      stop("candidate_loci must have 'chromosome' and 'position' columns")
    }
    
    if (verbose) message("  - Filtering for ", nrow(candidate_loci), " candidate loci")
    
    # Create temporary table for candidate loci
    temp_table_name <- paste0("temp_go_candidates_", as.integer(Sys.time()))
    
    DBI::dbExecute(con, paste0("CREATE TEMP TABLE ", temp_table_name, " (chromosome TEXT, position INTEGER)"))
    
    for (i in 1:nrow(candidate_loci)) {
      DBI::dbExecute(con, paste0("INSERT INTO ", temp_table_name, " VALUES (?, ?)"),
                    list(candidate_loci$chromosome[i], candidate_loci$position[i]))
    }
    
    loci_join <- paste0("JOIN flanking_sequences fs ON br.flanking_id = fs.flanking_id
                         JOIN vcf_data vd ON fs.vcf_id = vd.vcf_id
                         JOIN ", temp_table_name, " tc ON vd.chromosome = tc.chromosome AND vd.position = tc.position")
  } else {
    loci_join <- "JOIN flanking_sequences fs ON br.flanking_id = fs.flanking_id"
  }
  
  base_where <- if (length(base_conditions) > 0) paste("WHERE", paste(base_conditions, collapse = " AND ")) else ""
  
  # 1. Overview statistics by GO category
  if (verbose) message("  - Computing category overview...")
  
  overview_query <- paste0("
    SELECT 
      gt.go_category,
      COUNT(DISTINCT gt.go_id) as unique_go_terms,
      COUNT(gt.go_term_id) as total_annotations,
      COUNT(DISTINCT a.annotation_id) as annotated_sequences,
      COUNT(DISTINCT br.flanking_id) as annotated_loci
    FROM go_terms gt
    JOIN annotations a ON gt.annotation_id = a.annotation_id
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
    JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
    ", loci_join, "
    ", base_where, "
    GROUP BY gt.go_category
    ORDER BY gt.go_category
  ")
  
  if (length(params) > 0) {
    overview <- DBI::dbGetQuery(con, overview_query, params)
  } else {
    overview <- DBI::dbGetQuery(con, overview_query)
  }
  
  # 2. Top GO terms by category
  if (verbose) message("  - Finding top GO terms by category...")
  
  top_terms_query <- paste0("
    WITH term_counts AS (
      SELECT 
        gt.go_category,
        gt.go_id,
        gt.go_term,
        COUNT(gt.go_term_id) as frequency,
        COUNT(DISTINCT a.annotation_id) as unique_sequences
      FROM go_terms gt
      JOIN annotations a ON gt.annotation_id = a.annotation_id
      JOIN blast_results br ON a.blast_result_id = br.blast_result_id
      JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
      ", loci_join, "
      ", base_where, "
      GROUP BY gt.go_category, gt.go_id, gt.go_term
      HAVING COUNT(gt.go_term_id) >= ?
    ),
    ranked_terms AS (
      SELECT *,
        ROW_NUMBER() OVER (PARTITION BY go_category ORDER BY frequency DESC) as rank
      FROM term_counts
    )
    SELECT go_category, go_id, go_term, frequency, unique_sequences, rank
    FROM ranked_terms 
    WHERE rank <= 10
    ORDER BY go_category, rank
  ")
  
  top_terms_params <- c(params, list(min_frequency))
  top_terms <- DBI::dbGetQuery(con, top_terms_query, top_terms_params)
  
  # 3. Evidence code summary (if requested)
  evidence_summary <- NULL
  if (include_evidence) {
    if (verbose) message("  - Analyzing GO evidence codes...")
    
    evidence_query <- paste0("
      SELECT 
        gt.go_category,
        COALESCE(gt.go_evidence, 'Unknown') as evidence_code,
        COUNT(gt.go_term_id) as count,
        COUNT(DISTINCT gt.go_id) as unique_terms,
        ROUND(COUNT(gt.go_term_id) * 100.0 / SUM(COUNT(gt.go_term_id)) OVER (PARTITION BY gt.go_category), 2) as percent_in_category
      FROM go_terms gt
      JOIN annotations a ON gt.annotation_id = a.annotation_id
      JOIN blast_results br ON a.blast_result_id = br.blast_result_id
      JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
      ", loci_join, "
      ", base_where, "
      GROUP BY gt.go_category, gt.go_evidence
      ORDER BY gt.go_category, count DESC
    ")
    
    if (length(params) > 0) {
      evidence_summary <- DBI::dbGetQuery(con, evidence_query, params)
    } else {
      evidence_summary <- DBI::dbGetQuery(con, evidence_query)
    }
  }
  
  # 4. Coverage statistics across BLAST parameters
  if (verbose) message("  - Computing coverage statistics...")
  
  if (is.null(blast_param_id)) {
    coverage_query <- "
      SELECT 
        bp.blast_param_id,
        bp.blast_type,
        bp.db_name,
        COUNT(DISTINCT gt.go_id) as unique_go_terms,
        COUNT(gt.go_term_id) as total_annotations,
        COUNT(DISTINCT a.annotation_id) as annotated_sequences,
        COUNT(DISTINCT br.flanking_id) as annotated_loci
      FROM blast_parameters bp
      LEFT JOIN blast_results br ON bp.blast_param_id = br.blast_param_id
      LEFT JOIN annotations a ON br.blast_result_id = a.blast_result_id
      LEFT JOIN go_terms gt ON a.annotation_id = gt.annotation_id
      GROUP BY bp.blast_param_id, bp.blast_type, bp.db_name
      ORDER BY bp.blast_param_id
    "
    coverage_stats <- DBI::dbGetQuery(con, coverage_query)
  } else {
    # For specific blast_param_id, show comparison with total
    coverage_query <- "
      SELECT 
        'All BLAST parameters' as blast_type,
        COUNT(DISTINCT gt.go_id) as unique_go_terms,
        COUNT(gt.go_term_id) as total_annotations,
        COUNT(DISTINCT a.annotation_id) as annotated_sequences,
        COUNT(DISTINCT br.flanking_id) as annotated_loci
      FROM go_terms gt
      JOIN annotations a ON gt.annotation_id = a.annotation_id
      JOIN blast_results br ON a.blast_result_id = br.blast_result_id
      
      UNION ALL
      
      SELECT 
        'Selected BLAST parameter (' || bp.blast_param_id || ')' as blast_type,
        COUNT(DISTINCT gt.go_id) as unique_go_terms,
        COUNT(gt.go_term_id) as total_annotations,
        COUNT(DISTINCT a.annotation_id) as annotated_sequences,
        COUNT(DISTINCT br.flanking_id) as annotated_loci
      FROM go_terms gt
      JOIN annotations a ON gt.annotation_id = a.annotation_id
      JOIN blast_results br ON a.blast_result_id = br.blast_result_id
      JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
      WHERE bp.blast_param_id = ?
    "
    coverage_stats <- DBI::dbGetQuery(con, coverage_query, list(blast_param_id))
  }
  
  # 5. Functional module grouping (if requested)
  functional_modules <- NULL
  if (include_modules) {
    if (verbose) message("  - Grouping GO terms into functional modules...")
    functional_modules <- .create_go_functional_modules(top_terms, verbose = FALSE)
  }
  
  # Clean up temporary table if created
  if (!is.null(temp_table_name)) {
    DBI::dbExecute(con, paste0("DROP TABLE ", temp_table_name))
  }
  
  # Compile results
  result <- list(
    overview = overview,
    top_terms = top_terms,
    evidence_summary = evidence_summary,
    coverage_stats = coverage_stats,
    functional_modules = functional_modules,
    parameters = list(
      blast_param_id = blast_param_id,
      candidate_loci = !is.null(candidate_loci),
      include_evidence = include_evidence,
      include_modules = include_modules,
      min_frequency = min_frequency
    )
  )
  
  if (verbose) {
    message("GO annotation summary complete:")
    if (nrow(overview) > 0) {
      message("  - Found annotations in ", nrow(overview), " GO categories")
      message("  - Total unique GO terms: ", sum(overview$unique_go_terms))
      message("  - Total annotations: ", sum(overview$total_annotations))
    } else {
      message("  - No GO annotations found with current filters")
    }
  }
  
  return(result)
}

#' Summarize KEGG pathways in the database
#'
#' Provides comprehensive statistics about KEGG pathway annotations with optional
#' candidate loci filtering and functional module grouping by pathway categories.
#' Enhanced with KEGGREST and ggkegg integration for improved pathway analysis.
#'
#' @param con Database connection object
#' @param blast_param_id Integer. Specific BLAST parameter set to analyze. If NULL, analyzes all.
#' @param candidate_loci Data frame with chromosome and position columns, or NULL for all loci.
#' @param include_modules Logical. Group pathways into functional modules. Default is TRUE.
#' @param min_frequency Integer. Minimum frequency threshold for pathways to include. Default is 1.
#' @param use_keggrest_enhancement Logical. Enhance pathway names using KEGGREST. Default is FALSE.
#' @param use_ggkegg_analysis Logical. Include ggkegg-based network analysis. Default is FALSE.
#' @param verbose Logical. Print progress information. Default is TRUE.
#'
#' @return List containing KEGG pathway summaries:
#' \itemize{
#'   \item overview: Summary statistics of pathway coverage
#'   \item top_pathways: Most frequent pathways 
#'   \item pathway_modules: Pathways grouped by functional modules (if include_modules = TRUE)
#'   \item coverage_stats: Coverage statistics across BLAST parameters
#'   \item enhanced_analysis: KEGGREST/ggkegg analysis results (if requested)
#'   \item network_summary: Pathway network analysis (if use_ggkegg_analysis = TRUE)
#' }
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#' 
#' # Basic KEGG pathway analysis with module grouping
#' kegg_summary <- summarize_kegg_pathways(con, include_modules = TRUE)
#' print(kegg_summary$pathway_modules)
#' 
#' # Enhanced analysis with KEGGREST pathway name improvement
#' kegg_enhanced <- summarize_kegg_pathways(con, use_keggrest_enhancement = TRUE)
#' print(kegg_enhanced$enhanced_analysis)
#' 
#' # Full analysis with ggkegg network analysis for candidate loci
#' candidates <- data.frame(chromosome = c("LG1", "LG2"), position = c(12345, 67890))
#' kegg_full <- summarize_kegg_pathways(
#'   con, 
#'   candidate_loci = candidates,
#'   use_keggrest_enhancement = TRUE,
#'   use_ggkegg_analysis = TRUE
#' )
#' print(kegg_full$network_summary)
#' 
#' close_funseq_db(con)
#' }
#'
#' @export
summarize_kegg_pathways <- function(con, blast_param_id = NULL, candidate_loci = NULL,
                                   include_modules = TRUE, min_frequency = 1, 
                                   use_keggrest_enhancement = FALSE, use_ggkegg_analysis = FALSE,
                                   verbose = TRUE) {
  
  if (verbose) message("Summarizing KEGG pathways...")
  
  # Check if required tables exist
  tables <- DBI::dbListTables(con)
  required_tables <- c("kegg_references", "annotations", "blast_results", "blast_parameters")
  missing_tables <- required_tables[!required_tables %in% tables]
  
  if (length(missing_tables) > 0) {
    stop("Missing required tables: ", paste(missing_tables, collapse = ", "))
  }
  
  # Check if we have KEGG data
  kegg_count <- DBI::dbGetQuery(con, "SELECT COUNT(*) as count FROM kegg_references")$count
  if (kegg_count == 0) {
    warning("No KEGG pathway data found in database")
    return(list(
      overview = data.frame(),
      top_pathways = data.frame(),
      pathway_modules = NULL,
      coverage_stats = data.frame(),
      parameters = list(
        blast_param_id = blast_param_id,
        candidate_loci = !is.null(candidate_loci),
        include_modules = include_modules,
        min_frequency = min_frequency
      )
    ))
  }
  
  # Build base query with optional filters
  base_conditions <- c()
  params <- list()
  
  if (!is.null(blast_param_id)) {
    base_conditions <- c(base_conditions, "bp.blast_param_id = ?")
    params <- c(params, list(blast_param_id))
    if (verbose) message("  - Filtering for blast_param_id: ", blast_param_id)
  }
  
  # Handle candidate loci filtering
  temp_table_name <- NULL
  loci_join <- ""
  
  if (!is.null(candidate_loci)) {
    if (!all(c("chromosome", "position") %in% colnames(candidate_loci))) {
      stop("candidate_loci must have 'chromosome' and 'position' columns")
    }
    
    if (verbose) message("  - Filtering for ", nrow(candidate_loci), " candidate loci")
    
    # Create temporary table for candidate loci
    temp_table_name <- paste0("temp_kegg_candidates_", as.integer(Sys.time()))
    
    DBI::dbExecute(con, paste0("CREATE TEMP TABLE ", temp_table_name, " (chromosome TEXT, position INTEGER)"))
    
    for (i in 1:nrow(candidate_loci)) {
      DBI::dbExecute(con, paste0("INSERT INTO ", temp_table_name, " VALUES (?, ?)"),
                    list(candidate_loci$chromosome[i], candidate_loci$position[i]))
    }
    
    loci_join <- paste0("JOIN flanking_sequences fs ON br.flanking_id = fs.flanking_id
                         JOIN vcf_data vd ON fs.vcf_id = vd.vcf_id
                         JOIN ", temp_table_name, " tc ON vd.chromosome = tc.chromosome AND vd.position = tc.position")
  } else {
    loci_join <- "JOIN flanking_sequences fs ON br.flanking_id = fs.flanking_id"
  }
  
  base_where <- if (length(base_conditions) > 0) paste("WHERE", paste(base_conditions, collapse = " AND ")) else ""
  
  # 1. Overview statistics
  if (verbose) message("  - Computing pathway overview...")
  
  overview_query <- paste0("
    SELECT 
      COUNT(DISTINCT kr.kegg_id) as unique_pathways,
      COUNT(kr.kegg_ref_id) as total_pathway_assignments,
      COUNT(DISTINCT a.annotation_id) as annotated_sequences,
      COUNT(DISTINCT br.flanking_id) as annotated_loci
    FROM kegg_references kr
    JOIN annotations a ON kr.annotation_id = a.annotation_id
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
    JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
    ", loci_join, "
    ", base_where
  )
  
  if (length(params) > 0) {
    overview <- DBI::dbGetQuery(con, overview_query, params)
  } else {
    overview <- DBI::dbGetQuery(con, overview_query)
  }
  
  # 2. Top pathways by frequency
  if (verbose) message("  - Finding most frequent pathways...")
  
  top_pathways_query <- paste0("
    SELECT 
      kr.kegg_id,
      kr.pathway_name,
      COUNT(kr.kegg_ref_id) as frequency,
      COUNT(DISTINCT a.annotation_id) as unique_sequences,
      COUNT(DISTINCT br.flanking_id) as unique_loci
    FROM kegg_references kr
    JOIN annotations a ON kr.annotation_id = a.annotation_id
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
    JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
    ", loci_join, "
    ", base_where, "
    GROUP BY kr.kegg_id, kr.pathway_name
    HAVING COUNT(kr.kegg_ref_id) >= ?
    ORDER BY frequency DESC
    LIMIT 50
  ")
  
  top_pathways_params <- c(params, list(min_frequency))
  # Execute top pathways query with proper parameter handling
  num_placeholders <- nchar(top_pathways_query) - nchar(gsub("\\?", "", top_pathways_query))
  if (num_placeholders > 0 && length(top_pathways_params) == num_placeholders) {
    top_pathways <- DBI::dbGetQuery(con, top_pathways_query, top_pathways_params)
  } else if (num_placeholders == 0 && length(top_pathways_params) == 0) {
    top_pathways <- DBI::dbGetQuery(con, top_pathways_query)
  } else {
    stop("Parameter mismatch in top_pathways_query: ", num_placeholders, " placeholders, ", length(top_pathways_params), " parameters")
  }
  
  # 3. Coverage statistics across BLAST parameters
  if (verbose) message("  - Computing coverage statistics...")
  
  if (is.null(blast_param_id)) {
    coverage_query <- "
      SELECT 
        bp.blast_param_id,
        bp.blast_type,
        bp.db_name,
        COUNT(DISTINCT kr.kegg_id) as unique_pathways,
        COUNT(kr.kegg_ref_id) as total_assignments,
        COUNT(DISTINCT a.annotation_id) as annotated_sequences,
        COUNT(DISTINCT br.flanking_id) as annotated_loci
      FROM blast_parameters bp
      LEFT JOIN blast_results br ON bp.blast_param_id = br.blast_param_id
      LEFT JOIN annotations a ON br.blast_result_id = a.blast_result_id
      LEFT JOIN kegg_references kr ON a.annotation_id = kr.annotation_id
      GROUP BY bp.blast_param_id, bp.blast_type, bp.db_name
      ORDER BY bp.blast_param_id
    "
    coverage_stats <- DBI::dbGetQuery(con, coverage_query)
  } else {
    coverage_query <- "
      SELECT 
        'All BLAST parameters' as blast_type,
        COUNT(DISTINCT kr.kegg_id) as unique_pathways,
        COUNT(kr.kegg_ref_id) as total_assignments,
        COUNT(DISTINCT a.annotation_id) as annotated_sequences,
        COUNT(DISTINCT br.flanking_id) as annotated_loci
      FROM kegg_references kr
      JOIN annotations a ON kr.annotation_id = a.annotation_id
      JOIN blast_results br ON a.blast_result_id = br.blast_result_id
      
      UNION ALL
      
      SELECT 
        'Selected BLAST parameter (' || bp.blast_param_id || ')' as blast_type,
        COUNT(DISTINCT kr.kegg_id) as unique_pathways,
        COUNT(kr.kegg_ref_id) as total_assignments,
        COUNT(DISTINCT a.annotation_id) as annotated_sequences,
        COUNT(DISTINCT br.flanking_id) as annotated_loci
      FROM kegg_references kr
      JOIN annotations a ON kr.annotation_id = a.annotation_id
      JOIN blast_results br ON a.blast_result_id = br.blast_result_id
      JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
      WHERE bp.blast_param_id = ?
    "
    coverage_stats <- DBI::dbGetQuery(con, coverage_query, list(blast_param_id))
  }
  
  # 4. Pathway module grouping (if requested)
  pathway_modules <- NULL
  if (include_modules) {
    if (verbose) message("  - Grouping pathways into functional modules...")
    pathway_modules <- .create_kegg_functional_modules(top_pathways, verbose = FALSE)
  }
  
  # 5. Enhanced analysis using KEGGREST and ggkegg (if requested)
  enhanced_analysis <- NULL
  network_summary <- NULL
  
  if (use_keggrest_enhancement || use_ggkegg_analysis) {
    if (verbose) message("  - Performing enhanced KEGG analysis...")
    
    # KEGGREST enhancement for pathway names
    if (use_keggrest_enhancement) {
      if (verbose) message("    - Updating pathway names with KEGGREST...")
      
      # Only update if we have pathways with missing names
      missing_names_count <- sum(is.na(top_pathways$pathway_name) | 
                                top_pathways$pathway_name == "" | 
                                top_pathways$pathway_name == "-", na.rm = TRUE)
      
      if (missing_names_count > 0) {
        keggrest_results <- update_kegg_pathway_names_with_keggrest(
          con, 
          limit = 50, 
          verbose = verbose
        )
        enhanced_analysis$keggrest_update <- keggrest_results
      } else {
        if (verbose) message("      - No pathway names need updating")
        enhanced_analysis$keggrest_update <- list(message = "No pathway names needed updating")
      }
    }
    
    # ggkegg network analysis
    if (use_ggkegg_analysis) {
      if (verbose) message("    - Performing ggkegg pathway network analysis...")
      
      ggkegg_results <- analyze_kegg_modules(
        con, 
        candidate_loci = candidate_loci,
        blast_param_id = blast_param_id,
        include_network_analysis = TRUE,
        max_pathways_to_analyze = min(20, nrow(top_pathways)),
        verbose = verbose
      )
      
      network_summary <- list(
        pathway_networks = ggkegg_results$network_data,
        module_analysis = ggkegg_results$module_analysis,
        pathway_interactions = ggkegg_results$pathway_interactions,
        network_stats = ggkegg_results$coverage_stats
      )
      
      # Update pathway_modules with ggkegg results if available
      if (length(ggkegg_results$functional_modules) > 0) {
        enhanced_analysis$ggkegg_modules <- ggkegg_results$functional_modules
      }
    }
  }
  
  # Clean up temporary table if created
  if (!is.null(temp_table_name)) {
    DBI::dbExecute(con, paste0("DROP TABLE ", temp_table_name))
  }
  
  # Compile results
  result <- list(
    overview = overview,
    top_pathways = top_pathways,
    pathway_modules = pathway_modules,
    coverage_stats = coverage_stats,
    enhanced_analysis = enhanced_analysis,
    network_summary = network_summary,
    parameters = list(
      blast_param_id = blast_param_id,
      candidate_loci = !is.null(candidate_loci),
      include_modules = include_modules,
      min_frequency = min_frequency,
      use_keggrest_enhancement = use_keggrest_enhancement,
      use_ggkegg_analysis = use_ggkegg_analysis
    )
  )
  
  if (verbose) {
    message("KEGG pathway summary complete:")
    if (nrow(overview) > 0 && overview$total_pathway_assignments[1] > 0) {
      message("  - Unique pathways: ", overview$unique_pathways[1])
      message("  - Total assignments: ", overview$total_pathway_assignments[1])
      message("  - Annotated loci: ", overview$annotated_loci[1])
      
      # Enhanced analysis results
      if (use_keggrest_enhancement && !is.null(enhanced_analysis$keggrest_update)) {
        if (is.list(enhanced_analysis$keggrest_update) && "updated_pathways" %in% names(enhanced_analysis$keggrest_update)) {
          message("  - KEGGREST pathway names updated: ", enhanced_analysis$keggrest_update$updated_pathways)
        }
      }
      
      if (use_ggkegg_analysis && !is.null(network_summary)) {
        if (!is.null(network_summary$network_stats$pathways_with_networks)) {
          message("  - Pathway networks analyzed: ", network_summary$network_stats$pathways_with_networks)
        }
        if (!is.null(network_summary$pathway_interactions) && nrow(network_summary$pathway_interactions) > 0) {
          message("  - Pathway interactions found: ", nrow(network_summary$pathway_interactions))
        }
      }
    } else {
      message("  - No KEGG pathway assignments found with current filters")
    }
  }
  
  return(result)
}

#' Create comprehensive functional module summary
#'
#' Generates an integrated summary of functional modules combining GO terms, KEGG pathways,
#' and COG categories for a comprehensive functional characterization.
#'
#' @param con Database connection object
#' @param candidate_loci Data frame with chromosome and position columns, or NULL for all loci.
#' @param blast_param_id Integer. Specific BLAST parameter set to analyze. If NULL, analyzes all.
#' @param include_go Logical. Include GO term modules. Default is TRUE.
#' @param include_kegg Logical. Include KEGG pathway modules. Default is TRUE.
#' @param include_cog Logical. Include COG category analysis. Default is TRUE.
#' @param min_frequency Integer. Minimum frequency threshold for inclusion. Default is 2.
#' @param verbose Logical. Print progress information. Default is TRUE.
#'
#' @return List containing comprehensive functional module analysis:
#' \itemize{
#'   \item overview: Summary statistics across all annotation types
#'   \item go_modules: GO term functional modules (if include_go = TRUE)
#'   \item kegg_modules: KEGG pathway functional modules (if include_kegg = TRUE) 
#'   \item cog_modules: COG functional group analysis (if include_cog = TRUE)
#'   \item integrated_summary: Cross-platform functional summary
#'   \item module_overlap: Analysis of overlapping functional themes
#' }
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#' 
#' # Comprehensive module analysis for all loci
#' module_summary <- create_functional_module_summary(con)
#' print(module_summary$integrated_summary)
#' 
#' # Focus on candidate loci
#' candidates <- data.frame(chromosome = c("LG1", "LG2"), position = c(12345, 67890))
#' candidate_modules <- create_functional_module_summary(con, candidate_loci = candidates)
#' 
#' close_funseq_db(con)
#' }
#'
#' @export
create_functional_module_summary <- function(con, candidate_loci = NULL, blast_param_id = NULL,
                                           include_go = TRUE, include_kegg = TRUE, include_cog = TRUE,
                                           min_frequency = 2, verbose = TRUE) {
  
  if (verbose) message("Creating comprehensive functional module summary...")
  
  # Initialize results
  result <- list(
    overview = NULL,
    go_modules = NULL,
    kegg_modules = NULL,
    cog_modules = NULL,
    integrated_summary = NULL,
    module_overlap = NULL,
    parameters = list(
      candidate_loci = !is.null(candidate_loci),
      blast_param_id = blast_param_id,
      include_go = include_go,
      include_kegg = include_kegg,
      include_cog = include_cog,
      min_frequency = min_frequency
    )
  )
  
  # 1. GO modules analysis (if requested)
  if (include_go) {
    if (verbose) message("  - Analyzing GO functional modules...")
    go_summary <- summarize_go_annotations(
      con = con,
      blast_param_id = blast_param_id,
      candidate_loci = candidate_loci,
      include_modules = TRUE,
      min_frequency = min_frequency,
      verbose = FALSE
    )
    result$go_modules <- go_summary$functional_modules
  }
  
  # 2. KEGG modules analysis (if requested)
  if (include_kegg) {
    if (verbose) message("  - Analyzing KEGG pathway modules...")
    kegg_summary <- summarize_kegg_pathways(
      con = con,
      blast_param_id = blast_param_id,
      candidate_loci = candidate_loci,
      include_modules = TRUE,
      min_frequency = min_frequency,
      verbose = FALSE
    )
    result$kegg_modules <- kegg_summary$pathway_modules
  }
  
  # 3. COG modules analysis (if requested)
  if (include_cog) {
    if (verbose) message("  - Analyzing COG functional groups...")
    cog_summary <- summarize_cog_categories(
      con = con,
      blast_param_id = blast_param_id,
      include_functional_groups = TRUE,
      min_frequency = min_frequency,
      verbose = FALSE
    )
    result$cog_modules <- cog_summary$functional_groups
  }
  
  # 4. Create overview statistics
  if (verbose) message("  - Computing overview statistics...")
  
  overview_stats <- list()
  
  if (include_go && !is.null(result$go_modules) && length(result$go_modules) > 0) {
    overview_stats$go_module_count <- length(result$go_modules)
    overview_stats$go_total_terms <- sum(sapply(result$go_modules, function(m) m$term_count), na.rm = TRUE)
  }
  
  if (include_kegg && !is.null(result$kegg_modules) && length(result$kegg_modules) > 0) {
    overview_stats$kegg_module_count <- length(result$kegg_modules)
    overview_stats$kegg_total_pathways <- sum(sapply(result$kegg_modules, function(m) m$pathway_count), na.rm = TRUE)
  }
  
  if (include_cog && !is.null(result$cog_modules) && nrow(result$cog_modules) > 0) {
    overview_stats$cog_group_count <- nrow(result$cog_modules)
    overview_stats$cog_total_assignments <- sum(result$cog_modules$frequency, na.rm = TRUE)
  }
  
  result$overview <- overview_stats
  
  # 5. Create integrated summary
  if (verbose) message("  - Creating integrated functional summary...")
  result$integrated_summary <- .create_integrated_functional_summary(
    go_modules = result$go_modules,
    kegg_modules = result$kegg_modules,
    cog_modules = result$cog_modules,
    verbose = FALSE
  )
  
  # 6. Analyze module overlap
  if (verbose) message("  - Analyzing functional module overlap...")
  result$module_overlap <- .analyze_module_overlap(
    go_modules = result$go_modules,
    kegg_modules = result$kegg_modules,
    cog_modules = result$cog_modules,
    verbose = FALSE
  )
  
  if (verbose) {
    message("Functional module summary complete:")
    if (include_go && !is.null(result$go_modules)) {
      message("  - GO modules: ", length(result$go_modules))
    }
    if (include_kegg && !is.null(result$kegg_modules)) {
      message("  - KEGG modules: ", length(result$kegg_modules))
    }
    if (include_cog && !is.null(result$cog_modules)) {
      message("  - COG groups: ", nrow(result$cog_modules))
    }
    if (!is.null(result$integrated_summary)) {
      message("  - Integrated themes: ", length(result$integrated_summary))
    }
  }
  
  return(result)
}

#' Generate comprehensive candidate loci functional profile
#'
#' Creates a detailed functional characterization of candidate loci including
#' GO terms, KEGG pathways, protein annotations, and statistical comparisons
#' with background datasets.
#'
#' @param con Database connection object
#' @param candidate_loci Data frame with chromosome and position columns, or file_id of candidate VCF.
#' @param background_file_id Integer. File ID of background dataset for comparison. If NULL, auto-detects.
#' @param blast_param_id Integer. Specific BLAST parameter set to analyze. If NULL, analyzes all.
#' @param include_statistics Logical. Include statistical comparisons with background. Default is TRUE.
#' @param include_modules Logical. Include functional module analysis. Default is TRUE.
#' @param output_format Character. Format for results: "list", "tables", or "both". Default is "both".
#' @param verbose Logical. Print progress information. Default is TRUE.
#'
#' @return List containing comprehensive candidate loci profile:
#' \itemize{
#'   \item summary: Overview statistics and coverage metrics
#'   \item functional_modules: Integrated functional module analysis
#'   \item go_profile: GO term analysis specific to candidates
#'   \item kegg_profile: KEGG pathway analysis specific to candidates
#'   \item protein_profile: UniProt annotation analysis specific to candidates
#'   \item background_comparison: Statistical comparison with background (if include_statistics = TRUE)
#'   \item summary_tables: Publication-ready summary tables (if output_format includes "tables")
#' }
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#' 
#' # Profile specific candidate loci
#' candidates <- data.frame(
#'   chromosome = c("LG1", "LG2", "LG1"),
#'   position = c(12345, 67890, 23456)
#' )
#' profile <- generate_candidate_loci_profile(con, candidate_loci = candidates)
#' 
#' # Profile using candidate VCF file ID
#' profile <- generate_candidate_loci_profile(con, candidate_loci = 5)
#' 
#' print(profile$summary)
#' print(profile$functional_modules$integrated_summary)
#' 
#' close_funseq_db(con)
#' }
#'
#' @export
generate_candidate_loci_profile <- function(con, candidate_loci, background_file_id = NULL,
                                          blast_param_id = NULL, include_statistics = TRUE,
                                          include_modules = TRUE, output_format = "both",
                                          verbose = TRUE) {
  
  if (verbose) message("Generating comprehensive candidate loci profile...")
  
  # Validate output format
  if (!output_format %in% c("list", "tables", "both")) {
    stop("output_format must be 'list', 'tables', or 'both'")
  }
  
  # Handle different input types for candidate_loci
  if (is.numeric(candidate_loci) && length(candidate_loci) == 1) {
    # candidate_loci is a file_id
    candidate_file_id <- candidate_loci
    if (verbose) message("  - Using candidate file_id: ", candidate_file_id)
    
    # Extract loci coordinates from file
    loci_query <- "
      SELECT DISTINCT chromosome, position 
      FROM vcf_data 
      WHERE file_id = ?
      ORDER BY chromosome, position
    "
    candidate_coords <- DBI::dbGetQuery(con, loci_query, list(candidate_file_id))
    if (verbose) message("    - Extracted ", nrow(candidate_coords), " candidate coordinates")
    
  } else if (is.data.frame(candidate_loci)) {
    # candidate_loci is a data frame with coordinates
    if (!all(c("chromosome", "position") %in% colnames(candidate_loci))) {
      stop("candidate_loci data frame must have 'chromosome' and 'position' columns")
    }
    candidate_coords <- candidate_loci
    candidate_file_id <- NULL
    if (verbose) message("  - Using ", nrow(candidate_coords), " candidate loci coordinates")
    
  } else {
    stop("candidate_loci must be either a file_id (numeric) or data frame with chromosome/position columns")
  }
  
  # Auto-detect background file if not provided
  if (is.null(background_file_id)) {
    if (verbose) message("  - Auto-detecting background dataset...")
    
    input_files <- DBI::dbGetQuery(con, 
      "SELECT file_id, file_name FROM input_files WHERE file_type = 'vcf'")
    
    # Exclude candidate file if it's in the database
    if (!is.null(candidate_file_id)) {
      input_files <- input_files[input_files$file_id != candidate_file_id, ]
    }
    
    # Exclude files with "candidate" or similar in the name
    candidate_patterns <- c("candidate", "adaptive", "outlier", "fst", "selection")
    background_files <- input_files[!grepl(paste(candidate_patterns, collapse = "|"), 
                                          tolower(input_files$file_name)), ]
    
    if (nrow(background_files) == 0) {
      stop("No suitable background dataset found. Please specify background_file_id.")
    }
    
    background_file_id <- background_files$file_id[1]
    if (verbose) message("    - Using background file: ", background_files$file_name[1], " (ID: ", background_file_id, ")")
  }
  
  # Initialize results structure
  result <- list(
    summary = NULL,
    functional_modules = NULL,
    go_profile = NULL,
    kegg_profile = NULL,
    protein_profile = NULL,
    background_comparison = NULL,
    summary_tables = NULL,
    parameters = list(
      candidate_count = nrow(candidate_coords),
      background_file_id = background_file_id,
      blast_param_id = blast_param_id,
      include_statistics = include_statistics,
      include_modules = include_modules,
      output_format = output_format
    )
  )
  
  # 1. Get comprehensive functional profile for candidates
  if (verbose) message("  - Analyzing candidate functional profile...")
  
  candidate_profile <- get_functional_profile(
    con = con,
    candidate_loci = candidate_coords,
    blast_param_id = blast_param_id,
    include_blast_stats = TRUE,
    include_pathways = TRUE,
    verbose = TRUE
  )
  
  result$summary <- candidate_profile$loci_summary
  result$go_profile <- candidate_profile$go_profile
  result$protein_profile <- candidate_profile$protein_profile
  
  # Add KEGG profile if available
  if (!is.null(candidate_profile$pathway_profile)) {
    result$kegg_profile <- list(
      overview = data.frame(
        unique_pathways = length(unique(candidate_profile$pathway_profile$kegg_id)),
        total_assignments = nrow(candidate_profile$pathway_profile),
        annotated_loci = length(unique(candidate_profile$pathway_profile$unique_loci))
      ),
      top_pathways = candidate_profile$pathway_profile
    )
  }
  
  # 2. Create functional modules analysis (if requested)
  if (include_modules) {
    if (verbose) message("  - Creating functional module analysis...")
    
    result$functional_modules <- create_functional_module_summary(
      con = con,
      candidate_loci = candidate_coords,
      blast_param_id = blast_param_id,
      include_go = TRUE,
      include_kegg = TRUE,
      include_cog = TRUE,
      verbose = FALSE
    )
  }
  
  # 3. Statistical comparison with background (if requested)
  if (include_statistics) {
    if (verbose) message("  - Performing statistical comparison with background...")
    
    # Get background profile
    background_profile <- get_functional_profile(
      con = con,
      candidate_loci = NULL, # Use all loci for background
      blast_param_id = blast_param_id,
      include_blast_stats = TRUE,
      include_pathways = TRUE,
      verbose = FALSE
    )
    
    # Create comparison statistics
    comparison_stats <- list(
      coverage_comparison = list(
        candidate_coverage = if (result$summary$total_input_loci[1] > 0) {
          round(result$summary$annotated_loci[1] / result$summary$total_input_loci[1] * 100, 2)
        } else { 0 },
        background_coverage = if (background_profile$loci_summary$total_input_loci[1] > 0) {
          round(background_profile$loci_summary$annotated_loci[1] / background_profile$loci_summary$total_input_loci[1] * 100, 2)
        } else { 0 }
      ),
      annotation_enrichment = list(
        candidate_annotations_per_locus = if (result$summary$annotated_loci[1] > 0) {
          round(result$summary$total_annotations[1] / result$summary$annotated_loci[1], 2)
        } else { 0 },
        background_annotations_per_locus = if (background_profile$loci_summary$annotated_loci[1] > 0) {
          round(background_profile$loci_summary$total_annotations[1] / background_profile$loci_summary$annotated_loci[1], 2)
        } else { 0 }
      ),
      functional_diversity = list(
        candidate_go_diversity = if (!is.null(result$go_profile$overview)) {
          result$go_profile$overview$unique_go_terms[1]
        } else { 0 },
        background_go_diversity = if (!is.null(background_profile$go_profile$overview)) {
          background_profile$go_profile$overview$unique_go_terms[1]
        } else { 0 }
      )
    )
    
    result$background_comparison <- comparison_stats
  }
  
  # 4. Create summary tables (if requested)
  if (output_format %in% c("tables", "both")) {
    if (verbose) message("  - Creating summary tables...")
    
    summary_tables <- list()
    
    # Overview table
    overview_table <- data.frame(
      Metric = c(
        "Total Candidate Loci",
        "Annotated Loci",
        "Annotation Coverage (%)",
        "Total Functional Annotations",
        "Unique Proteins",
        "Unique GO Terms",
        "Unique KEGG Pathways"
      ),
      Value = c(
        result$summary$total_input_loci[1],
        result$summary$annotated_loci[1],
        result$summary$annotation_coverage_percent[1],
        result$summary$total_annotations[1],
        result$summary$unique_proteins[1],
        if (!is.null(result$go_profile$overview)) result$go_profile$overview$unique_go_terms[1] else 0,
        if (!is.null(result$kegg_profile)) result$kegg_profile$overview$unique_pathways[1] else 0
      )
    )
    
    summary_tables$overview <- overview_table
    
    # Top GO terms table (if available)
    if (!is.null(result$go_profile$top_terms) && nrow(result$go_profile$top_terms) > 0) {
      top_go_table <- result$go_profile$top_terms[1:min(10, nrow(result$go_profile$top_terms)), 
                                                  c("go_id", "go_term", "go_category", "frequency")]
      colnames(top_go_table) <- c("GO ID", "GO Term", "Category", "Frequency")
      summary_tables$top_go_terms <- top_go_table
    }
    
    # Top KEGG pathways table (if available)
    if (!is.null(result$kegg_profile$top_pathways) && nrow(result$kegg_profile$top_pathways) > 0) {
      top_kegg_table <- result$kegg_profile$top_pathways[1:min(10, nrow(result$kegg_profile$top_pathways)), 
                                                         c("kegg_id", "pathway_name", "frequency")]
      colnames(top_kegg_table) <- c("KEGG ID", "Pathway Name", "Frequency")
      summary_tables$top_kegg_pathways <- top_kegg_table
    }
    
    # Functional modules summary (if available)
    if (!is.null(result$functional_modules$integrated_summary)) {
      modules_summary <- data.frame(
        Functional_Theme = names(result$functional_modules$integrated_summary),
        Sources = sapply(result$functional_modules$integrated_summary, function(theme) {
          paste(names(theme$sources), collapse = ", ")
        }),
        Description = sapply(result$functional_modules$integrated_summary, function(theme) {
          theme$description
        })
      )
      summary_tables$functional_modules <- modules_summary
    }
    
    result$summary_tables <- summary_tables
  }
  
  if (verbose) {
    message("Candidate loci profile complete:")
    message("  - Input loci: ", result$summary$total_input_loci[1])
    message("  - Annotated loci: ", result$summary$annotated_loci[1], 
            " (", result$summary$annotation_coverage_percent[1], "%)")
    message("  - Functional annotations: ", result$summary$total_annotations[1])
    if (include_modules && !is.null(result$functional_modules)) {
      message("  - Functional themes identified: ", length(result$functional_modules$integrated_summary))
    }
  }
  
  return(result)
}

#' Create functional comparison tables between datasets
#'
#' Generates publication-ready comparison tables showing functional differences
#' between multiple datasets (e.g., candidate vs background, or multiple candidate sets).
#'
#' @param con Database connection object
#' @param dataset_list List. Named list of datasets, where each element is either:
#'   - Data frame with chromosome/position columns
#'   - Numeric file_id
#'   - NULL for all loci
#' @param blast_param_id Integer. Specific BLAST parameter set to analyze. If NULL, analyzes all.
#' @param comparison_metrics Character vector. Metrics to compare: c("coverage", "diversity", "modules", "top_terms").
#' @param max_terms Integer. Maximum number of top terms to show per dataset. Default is 10.
#' @param include_statistics Logical. Include statistical significance tests. Default is TRUE.
#' @param verbose Logical. Print progress information. Default is TRUE.
#'
#' @return List containing comparison tables:
#' \itemize{
#'   \item overview_comparison: Side-by-side overview statistics
#'   \item coverage_comparison: Annotation coverage comparison
#'   \item diversity_comparison: Functional diversity metrics comparison
#'   \item top_terms_comparison: Most frequent terms per dataset
#'   \item module_comparison: Functional module representation comparison (if "modules" in comparison_metrics)
#'   \item statistical_tests: Statistical significance tests (if include_statistics = TRUE)
#' }
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#' 
#' # Compare candidate vs background
#' datasets <- list(
#'   "Candidates" = data.frame(chromosome = c("LG1", "LG2"), position = c(12345, 67890)),
#'   "Background" = NULL  # All loci
#' )
#' comparison <- create_functional_comparison_tables(con, datasets)
#' print(comparison$overview_comparison)
#' 
#' # Compare multiple candidate sets
#' datasets <- list(
#'   "Adaptive_Set1" = 5,   # file_id
#'   "Adaptive_Set2" = 6,   # file_id
#'   "Neutral" = NULL       # All loci
#' )
#' comparison <- create_functional_comparison_tables(con, datasets)
#' 
#' close_funseq_db(con)
#' }
#'
#' @export
create_functional_comparison_tables <- function(con, dataset_list, blast_param_id = NULL,
                                              comparison_metrics = c("coverage", "diversity", "top_terms"),
                                              max_terms = 10, include_statistics = TRUE, verbose = TRUE) {
  
  if (verbose) message("Creating functional comparison tables...")
  
  # Validate inputs
  if (!is.list(dataset_list) || is.null(names(dataset_list))) {
    stop("dataset_list must be a named list")
  }
  
  if (length(dataset_list) < 2) {
    stop("dataset_list must contain at least 2 datasets for comparison")
  }
  
  valid_metrics <- c("coverage", "diversity", "modules", "top_terms")
  if (!all(comparison_metrics %in% valid_metrics)) {
    stop("comparison_metrics must be from: ", paste(valid_metrics, collapse = ", "))
  }
  
  # Initialize results
  result <- list(
    overview_comparison = NULL,
    coverage_comparison = NULL,
    diversity_comparison = NULL,
    top_terms_comparison = NULL,
    module_comparison = NULL,
    statistical_tests = NULL,
    parameters = list(
      dataset_names = names(dataset_list),
      dataset_count = length(dataset_list),
      blast_param_id = blast_param_id,
      comparison_metrics = comparison_metrics,
      max_terms = max_terms,
      include_statistics = include_statistics
    )
  )
  
  # Process each dataset
  if (verbose) message("  - Processing ", length(dataset_list), " datasets...")
  
  dataset_profiles <- list()
  
  for (dataset_name in names(dataset_list)) {
    if (verbose) message("    - Processing: ", dataset_name)
    
    dataset_data <- dataset_list[[dataset_name]]
    
    # Convert dataset specification to candidate_loci format
    if (is.null(dataset_data)) {
      # All loci
      candidate_loci <- NULL
    } else if (is.numeric(dataset_data) && length(dataset_data) == 1) {
      # File ID - extract coordinates
      loci_query <- "SELECT DISTINCT chromosome, position FROM vcf_data WHERE file_id = ?"
      candidate_loci <- DBI::dbGetQuery(con, loci_query, list(dataset_data))
    } else if (is.data.frame(dataset_data)) {
      # Coordinate data frame
      if (!all(c("chromosome", "position") %in% colnames(dataset_data))) {
        stop("Dataset '", dataset_name, "' must have chromosome and position columns")
      }
      candidate_loci <- dataset_data
    } else {
      stop("Invalid dataset specification for '", dataset_name, "'")
    }
    
    # Get functional profile
    profile <- get_functional_profile(
      con = con,
      candidate_loci = candidate_loci,
      blast_param_id = blast_param_id,
      include_blast_stats = TRUE,
      include_pathways = TRUE,
      verbose = FALSE
    )
    
    dataset_profiles[[dataset_name]] <- profile
  }
  
  # 1. Overview comparison
  if ("coverage" %in% comparison_metrics || "diversity" %in% comparison_metrics) {
    if (verbose) message("  - Creating overview comparison...")
    
    overview_data <- data.frame(
      Dataset = names(dataset_profiles),
      Total_Loci = sapply(dataset_profiles, function(p) p$loci_summary$total_input_loci[1]),
      Annotated_Loci = sapply(dataset_profiles, function(p) p$loci_summary$annotated_loci[1]),
      Coverage_Percent = sapply(dataset_profiles, function(p) p$loci_summary$annotation_coverage_percent[1]),
      Total_Annotations = sapply(dataset_profiles, function(p) p$loci_summary$total_annotations[1]),
      Unique_Proteins = sapply(dataset_profiles, function(p) p$loci_summary$unique_proteins[1]),
      Unique_GO_Terms = sapply(dataset_profiles, function(p) {
        if (!is.null(p$go_profile$overview)) p$go_profile$overview$unique_go_terms[1] else 0
      }),
      stringsAsFactors = FALSE
    )
    
    result$overview_comparison <- overview_data
  }
  
  # 2. Coverage comparison
  if ("coverage" %in% comparison_metrics) {
    if (verbose) message("  - Creating coverage comparison...")
    
    coverage_data <- data.frame(
      Dataset = names(dataset_profiles),
      Annotation_Coverage = sapply(dataset_profiles, function(p) p$loci_summary$annotation_coverage_percent[1]),
      Annotations_Per_Locus = sapply(dataset_profiles, function(p) {
        if (p$loci_summary$annotated_loci[1] > 0) {
          round(p$loci_summary$total_annotations[1] / p$loci_summary$annotated_loci[1], 2)
        } else { 0 }
      }),
      GO_Coverage = sapply(dataset_profiles, function(p) {
        if (!is.null(p$go_profile$coverage_stats)) {
          p$go_profile$coverage_stats$go_coverage_percent[1]
        } else { 0 }
      }),
      stringsAsFactors = FALSE
    )
    
    result$coverage_comparison <- coverage_data
  }
  
  # 3. Diversity comparison
  if ("diversity" %in% comparison_metrics) {
    if (verbose) message("  - Creating diversity comparison...")
    
    diversity_data <- data.frame(
      Dataset = names(dataset_profiles),
      Functional_Diversity = sapply(dataset_profiles, function(p) {
        if (!is.null(p$functional_diversity)) {
          p$functional_diversity$go_terms_per_locus[1]
        } else { 0 }
      }),
      Protein_Diversity = sapply(dataset_profiles, function(p) {
        if (!is.null(p$functional_diversity)) {
          p$functional_diversity$proteins_per_locus[1]
        } else { 0 }
      }),
      GO_Categories = sapply(dataset_profiles, function(p) {
        if (!is.null(p$functional_diversity)) {
          p$functional_diversity$go_categories[1]
        } else { 0 }
      }),
      stringsAsFactors = FALSE
    )
    
    result$diversity_comparison <- diversity_data
  }
  
  # 4. Top terms comparison
  if ("top_terms" %in% comparison_metrics) {
    if (verbose) message("  - Creating top terms comparison...")
    
    # Combine top GO terms from all datasets
    all_go_terms <- list()
    
    for (dataset_name in names(dataset_profiles)) {
      profile <- dataset_profiles[[dataset_name]]
      if (!is.null(profile$go_profile$top_terms) && nrow(profile$go_profile$top_terms) > 0) {
        top_terms <- profile$go_profile$top_terms[1:min(max_terms, nrow(profile$go_profile$top_terms)), ]
        top_terms$Dataset <- dataset_name
        all_go_terms[[dataset_name]] <- top_terms[, c("Dataset", "go_id", "go_term", "go_category", "frequency")]
      }
    }
    
    if (length(all_go_terms) > 0) {
      result$top_terms_comparison <- do.call(rbind, all_go_terms)
      rownames(result$top_terms_comparison) <- NULL
    }
  }
  
  # 5. Module comparison (if requested)
  if ("modules" %in% comparison_metrics) {
    if (verbose) message("  - Creating functional module comparison...")
    
    module_summaries <- list()
    
    for (dataset_name in names(dataset_profiles)) {
      profile <- dataset_profiles[[dataset_name]]
      
      # Get candidate loci for this dataset
      dataset_data <- dataset_list[[dataset_name]]
      if (is.null(dataset_data)) {
        candidate_loci <- NULL
      } else if (is.numeric(dataset_data) && length(dataset_data) == 1) {
        loci_query <- "SELECT DISTINCT chromosome, position FROM vcf_data WHERE file_id = ?"
        candidate_loci <- DBI::dbGetQuery(con, loci_query, list(dataset_data))
      } else {
        candidate_loci <- dataset_data
      }
      
      # Get module analysis
      module_analysis <- create_functional_module_summary(
        con = con,
        candidate_loci = candidate_loci,
        blast_param_id = blast_param_id,
        verbose = FALSE
      )
      
      if (!is.null(module_analysis$integrated_summary)) {
        module_summary <- data.frame(
          Dataset = dataset_name,
          Functional_Theme = names(module_analysis$integrated_summary),
          Sources = sapply(module_analysis$integrated_summary, function(theme) {
            paste(names(theme$sources), collapse = ", ")
          }),
          stringsAsFactors = FALSE
        )
        module_summaries[[dataset_name]] <- module_summary
      }
    }
    
    if (length(module_summaries) > 0) {
      result$module_comparison <- do.call(rbind, module_summaries)
      rownames(result$module_comparison) <- NULL
    }
  }
  
  # 6. Statistical tests (if requested)
  if (include_statistics && length(dataset_profiles) == 2) {
    if (verbose) message("  - Performing statistical tests...")
    
    dataset_names <- names(dataset_profiles)
    profile1 <- dataset_profiles[[dataset_names[1]]]
    profile2 <- dataset_profiles[[dataset_names[2]]]
    
    # Coverage comparison test
    coverage_test <- NULL
    if (profile1$loci_summary$total_input_loci[1] > 0 && profile2$loci_summary$total_input_loci[1] > 0) {
      # Proportion test for annotation coverage
      coverage_test <- prop.test(
        x = c(profile1$loci_summary$annotated_loci[1], profile2$loci_summary$annotated_loci[1]),
        n = c(profile1$loci_summary$total_input_loci[1], profile2$loci_summary$total_input_loci[1])
      )
    }
    
    # Annotation density comparison (if both have annotated loci)
    annotation_test <- NULL
    if (profile1$loci_summary$annotated_loci[1] > 0 && profile2$loci_summary$annotated_loci[1] > 0) {
      # t-test for annotations per locus (requires individual locus data)
      # For now, just report the comparison metrics
      annotation_test <- list(
        annotations_per_locus_1 = round(profile1$loci_summary$total_annotations[1] / profile1$loci_summary$annotated_loci[1], 2),
        annotations_per_locus_2 = round(profile2$loci_summary$total_annotations[1] / profile2$loci_summary$annotated_loci[1], 2)
      )
    }
    
    result$statistical_tests <- list(
      coverage_test = coverage_test,
      annotation_test = annotation_test,
      comparison_names = dataset_names
    )
  }
  
  if (verbose) {
    message("Functional comparison complete:")
    message("  - Datasets compared: ", length(dataset_profiles))
    message("  - Comparison metrics: ", paste(comparison_metrics, collapse = ", "))
    if (include_statistics) {
      message("  - Statistical tests performed: ", !is.null(result$statistical_tests))
    }
  }
  
  return(result)
}

#' Summarize UniProt annotations in the database
#'
#' Provides comprehensive statistics about UniProt protein annotations
#' including gene names, protein functions, and annotation coverage.
#'
#' @param con Database connection object
#' @param blast_param_id Integer. Specific BLAST parameter set to analyze. If NULL, analyzes all.
#' @param include_gene_names Logical. Include gene name analysis. Default is TRUE.
#' @param top_n Integer. Number of top items to include in summaries. Default is 20.
#' @param verbose Logical. Print progress information. Default is TRUE.
#'
#' @return List containing UniProt annotation summaries:
#' \itemize{
#'   \item overview: Basic statistics about UniProt annotations
#'   \item top_proteins: Most frequently annotated proteins
#'   \item gene_name_stats: Gene name coverage and distribution (if include_gene_names = TRUE)
#'   \item annotation_quality: Quality metrics for annotations
#' }
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#' 
#' # Analyze all UniProt annotations
#' uniprot_summary <- summarize_uniprot_annotations(con)
#' print(uniprot_summary$overview)
#' 
#' # Focus on specific BLAST results
#' uniprot_summary_blast1 <- summarize_uniprot_annotations(con, blast_param_id = 1)
#' 
#' close_funseq_db(con)
#' }
#'
#' @export
summarize_uniprot_annotations <- function(con, blast_param_id = NULL, include_gene_names = TRUE, 
                                        top_n = 20, verbose = TRUE) {
  
  if (verbose) message("Summarizing UniProt annotations...")
  
  # Check if required tables exist
  tables <- DBI::dbListTables(con)
  required_tables <- c("annotations", "blast_results", "blast_parameters")
  missing_tables <- required_tables[!required_tables %in% tables]
  
  if (length(missing_tables) > 0) {
    stop("Missing required tables: ", paste(missing_tables, collapse = ", "))
  }
  
  # Build base query with optional blast_param_id filter
  base_where <- ""
  params <- list()
  
  if (!is.null(blast_param_id)) {
    base_where <- "WHERE bp.blast_param_id = ?"
    params <- list(blast_param_id)
    if (verbose) message("  - Filtering for blast_param_id: ", blast_param_id)
  }
  
  # 1. Overview statistics
  if (verbose) message("  - Computing overview statistics...")
  
  overview_query <- paste0("
    SELECT 
      COUNT(DISTINCT a.uniprot_accession) as unique_uniprot_entries,
      COUNT(a.annotation_id) as total_annotations,
      COUNT(DISTINCT br.flanking_id) as annotated_loci,
      COUNT(CASE WHEN a.entry_name IS NOT NULL AND a.entry_name != '' THEN 1 END) as entries_with_names,
      COUNT(CASE WHEN a.gene_names IS NOT NULL AND a.gene_names != '' THEN 1 END) as entries_with_gene_names,
      ROUND(COUNT(CASE WHEN a.entry_name IS NOT NULL AND a.entry_name != '' THEN 1 END) * 100.0 / COUNT(a.annotation_id), 2) as percent_with_names,
      ROUND(COUNT(CASE WHEN a.gene_names IS NOT NULL AND a.gene_names != '' THEN 1 END) * 100.0 / COUNT(a.annotation_id), 2) as percent_with_gene_names
    FROM annotations a
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
    JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
    ", base_where
  )
  
  # Execute overview query with proper parameter handling
  num_placeholders <- nchar(overview_query) - nchar(gsub("\\?", "", overview_query))
  if (num_placeholders > 0 && length(params) == num_placeholders) {
    overview <- DBI::dbGetQuery(con, overview_query, params)
  } else if (num_placeholders == 0 && length(params) == 0) {
    overview <- DBI::dbGetQuery(con, overview_query)
  } else {
    stop("Parameter mismatch in overview_query: ", num_placeholders, " placeholders, ", length(params), " parameters")
  }
  
  # 2. Top proteins by frequency
  if (verbose) message("  - Finding most frequent protein annotations...")
  
  top_proteins_query <- paste0("
    SELECT 
      a.uniprot_accession,
      a.entry_name,
      a.gene_names,
      COUNT(a.annotation_id) as frequency,
      COUNT(DISTINCT br.flanking_id) as unique_loci
    FROM annotations a
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
    JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
    ", base_where, "
    GROUP BY a.uniprot_accession, a.entry_name, a.gene_names
    ORDER BY frequency DESC
    LIMIT ?
  ")
  
  top_proteins_params <- c(params, list(top_n))
  # Execute top proteins query with proper parameter handling
  num_placeholders <- nchar(top_proteins_query) - nchar(gsub("\\?", "", top_proteins_query))
  if (num_placeholders > 0 && length(top_proteins_params) == num_placeholders) {
    top_proteins <- DBI::dbGetQuery(con, top_proteins_query, top_proteins_params)
  } else {
    stop("Parameter mismatch in top_proteins_query: ", num_placeholders, " placeholders, ", length(top_proteins_params), " parameters")
  }
  
  # 3. Gene name statistics (if requested)
  gene_name_stats <- NULL
  if (include_gene_names) {
    if (verbose) message("  - Analyzing gene name patterns...")
    
    gene_stats_query <- paste0("
      WITH gene_names_split AS (
        SELECT 
          a.annotation_id,
          TRIM(SUBSTR(a.gene_names, 1, CASE 
            WHEN INSTR(a.gene_names, ';') > 0 THEN INSTR(a.gene_names, ';') - 1 
            ELSE LENGTH(a.gene_names) 
          END)) as primary_gene_name,
          a.gene_names as full_gene_names
        FROM annotations a
        JOIN blast_results br ON a.blast_result_id = br.blast_result_id
        JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
        ", base_where, "
        AND a.gene_names IS NOT NULL AND a.gene_names != ''
      ),
      gene_name_counts AS (
        SELECT 
          primary_gene_name,
          COUNT(*) as frequency,
          COUNT(DISTINCT annotation_id) as unique_annotations
        FROM gene_names_split
        WHERE primary_gene_name != ''
        GROUP BY primary_gene_name
      )
      SELECT 
        primary_gene_name,
        frequency,
        unique_annotations,
        ROUND(frequency * 100.0 / SUM(frequency) OVER (), 2) as percent_of_total
      FROM gene_name_counts
      ORDER BY frequency DESC
      LIMIT ?
    ")
    
    gene_stats_params <- c(params, list(top_n))
    # Execute gene stats query with proper parameter handling
    num_placeholders <- nchar(gene_stats_query) - nchar(gsub("\\?", "", gene_stats_query))
    if (num_placeholders > 0 && length(gene_stats_params) == num_placeholders) {
      gene_name_stats <- DBI::dbGetQuery(con, gene_stats_query, gene_stats_params)
    } else if (num_placeholders == 0 && length(gene_stats_params) == 0) {
      gene_name_stats <- DBI::dbGetQuery(con, gene_stats_query)
    } else {
      stop("Parameter mismatch in gene_stats_query: ", num_placeholders, " placeholders, ", length(gene_stats_params), " parameters")
    }
  }
  
  # 4. Annotation quality metrics
  if (verbose) message("  - Computing annotation quality metrics...")
  
  quality_query <- paste0("
    SELECT 
      'Total annotations' as metric,
      COUNT(a.annotation_id) as count,
      '' as description
    FROM annotations a
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
    JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
    ", base_where, "
    
    UNION ALL
    
    SELECT 
      'Complete annotations' as metric,
      COUNT(a.annotation_id) as count,
      'Have both entry_name and gene_names' as description
    FROM annotations a
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
    JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
    ", base_where, "
    AND a.entry_name IS NOT NULL AND a.entry_name != ''
    AND a.gene_names IS NOT NULL AND a.gene_names != ''
    
    UNION ALL
    
    SELECT 
      'Partial annotations' as metric,
      COUNT(a.annotation_id) as count,
      'Have either entry_name OR gene_names' as description
    FROM annotations a
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
    JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
    ", base_where, "
    AND ((a.entry_name IS NOT NULL AND a.entry_name != '') OR 
         (a.gene_names IS NOT NULL AND a.gene_names != ''))
    AND NOT (a.entry_name IS NOT NULL AND a.entry_name != '' AND 
             a.gene_names IS NOT NULL AND a.gene_names != '')
    
    UNION ALL
    
    SELECT 
      'Minimal annotations' as metric,
      COUNT(a.annotation_id) as count,
      'Only UniProt accession available' as description
    FROM annotations a
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
    JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
    ", base_where, "
    AND (a.entry_name IS NULL OR a.entry_name = '')
    AND (a.gene_names IS NULL OR a.gene_names = '')
  ")
  
  # Execute quality query with proper parameter handling
  num_placeholders <- nchar(quality_query) - nchar(gsub("\\?", "", quality_query))
  if (num_placeholders > 0 && length(params) == num_placeholders) {
    annotation_quality <- DBI::dbGetQuery(con, quality_query, params)
  } else if (num_placeholders == 0 && length(params) == 0) {
    annotation_quality <- DBI::dbGetQuery(con, quality_query)
  } else {
    stop("Parameter mismatch in quality_query: ", num_placeholders, " placeholders, ", length(params), " parameters")
  }
  
  # Compile results
  result <- list(
    overview = overview,
    top_proteins = top_proteins,
    gene_name_stats = gene_name_stats,
    annotation_quality = annotation_quality,
    parameters = list(
      blast_param_id = blast_param_id,
      include_gene_names = include_gene_names,
      top_n = top_n
    )
  )
  
  if (verbose) {
    message("UniProt annotation summary complete:")
    if (nrow(overview) > 0 && overview$total_annotations[1] > 0) {
      message("  - Unique UniProt entries: ", overview$unique_uniprot_entries[1])
      message("  - Total annotations: ", overview$total_annotations[1])
      message("  - Annotated loci: ", overview$annotated_loci[1])
      message("  - Coverage: ", overview$percent_with_names[1], "% have entry names, ",
              overview$percent_with_gene_names[1], "% have gene names")
    } else {
      message("  - No UniProt annotations found with current filters")
    }
  }
  
  return(result)
}

#' Get comprehensive functional profile of candidate loci
#'
#' Provides a complete functional characterization combining GO terms, UniProt annotations,
#' KEGG pathways, and BLAST hit statistics for a set of candidate loci.
#'
#' @param con Database connection object
#' @param candidate_loci Data frame with chromosome and position columns, or NULL for all annotated loci
#' @param blast_param_id Integer. Specific BLAST parameter set to use. If NULL, uses all available.
#' @param include_blast_stats Logical. Include BLAST hit quality statistics. Default is TRUE.
#' @param include_pathways Logical. Include KEGG pathway information. Default is TRUE.
#' @param verbose Logical. Print progress information. Default is TRUE.
#'
#' @return List containing comprehensive functional profile:
#' \itemize{
#'   \item loci_summary: Overview of input loci and annotation coverage
#'   \item go_profile: GO term distribution and enriched categories
#'   \item protein_profile: UniProt annotation summary
#'   \item pathway_profile: KEGG pathway annotations (if include_pathways = TRUE)
#'   \item blast_quality: BLAST hit statistics (if include_blast_stats = TRUE)
#'   \item functional_diversity: Measures of functional diversity
#' }
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#' 
#' # Profile all annotated loci
#' full_profile <- get_functional_profile(con)
#' 
#' # Profile specific candidate loci
#' candidates <- data.frame(
#'   chromosome = c("LG1", "LG2", "LG1"),
#'   position = c(12345, 67890, 23456)
#' )
#' candidate_profile <- get_functional_profile(con, candidate_loci = candidates)
#' 
#' close_funseq_db(con)
#' }
#'
#' @export
get_functional_profile <- function(con, candidate_loci = NULL, blast_param_id = NULL, 
                                 include_blast_stats = TRUE, include_pathways = TRUE, 
                                 verbose = TRUE) {
  
  if (verbose) message("Generating comprehensive functional profile...")
  
  # Check if required tables exist
  tables <- DBI::dbListTables(con)
  required_tables <- c("vcf_data", "flanking_sequences", "blast_results", "annotations")
  missing_tables <- required_tables[!required_tables %in% tables]
  
  if (length(missing_tables) > 0) {
    stop("Missing required tables: ", paste(missing_tables, collapse = ", "))
  }
  
  # Build base filtering conditions
  base_conditions <- c()
  params <- list()
  
  if (!is.null(blast_param_id)) {
    base_conditions <- c(base_conditions, "bp.blast_param_id = ?")
    params <- c(params, list(blast_param_id))
    if (verbose) message("  - Using blast_param_id: ", blast_param_id)
  } else {
    if (verbose) message("  - No blast_param_id filter applied")
  }
  
  # Handle candidate loci filtering
  loci_condition <- ""
  if (!is.null(candidate_loci)) {
    if (!all(c("chromosome", "position") %in% colnames(candidate_loci))) {
      stop("candidate_loci must have 'chromosome' and 'position' columns")
    }
    
    if (verbose) message("  - Filtering for ", nrow(candidate_loci), " candidate loci")
    
    # Create temporary table for candidate loci
    temp_table_name <- paste0("temp_candidates_", as.integer(Sys.time()))
    
    # Insert candidate loci into temporary table
    DBI::dbExecute(con, paste0("CREATE TEMP TABLE ", temp_table_name, " (chromosome TEXT, position INTEGER)"))
    
    for (i in 1:nrow(candidate_loci)) {
      DBI::dbExecute(con, paste0("INSERT INTO ", temp_table_name, " VALUES (?, ?)"),
                    list(candidate_loci$chromosome[i], candidate_loci$position[i]))
    }
    
    loci_condition <- paste0("AND EXISTS (SELECT 1 FROM ", temp_table_name, " tc WHERE vd.chromosome = tc.chromosome AND vd.position = tc.position)")
  }
  
  # Combine all conditions
  where_clause <- ""
  all_conditions <- c()
  
  if (length(base_conditions) > 0) {
    all_conditions <- c(all_conditions, base_conditions)
  }
  
  if (loci_condition != "") {
    # Remove "AND " from the beginning of loci_condition
    loci_condition_clean <- substr(loci_condition, 5, nchar(loci_condition))
    all_conditions <- c(all_conditions, loci_condition_clean)
  }
  
  if (length(all_conditions) > 0) {
    where_clause <- paste("WHERE", paste(all_conditions, collapse = " AND "))
  }
  
  if (verbose) message("  - WHERE clause: ", if (where_clause == "") "NONE" else where_clause)
  if (verbose) message("  - Final params list length: ", length(params))
  
  # 1. Loci summary
  if (verbose) message("  - Computing loci summary...")
  
  loci_summary_query <- paste0("
    SELECT 
      COUNT(DISTINCT vd.vcf_id) as total_input_loci,
      COUNT(DISTINCT CASE WHEN a.annotation_id IS NOT NULL THEN vd.vcf_id END) as annotated_loci,
      COUNT(DISTINCT a.annotation_id) as total_annotations,
      COUNT(DISTINCT a.uniprot_accession) as unique_proteins,
      ROUND(COUNT(DISTINCT CASE WHEN a.annotation_id IS NOT NULL THEN vd.vcf_id END) * 100.0 / 
            COUNT(DISTINCT vd.vcf_id), 2) as annotation_coverage_percent
    FROM vcf_data vd
    LEFT JOIN flanking_sequences fs ON vd.vcf_id = fs.vcf_id
    LEFT JOIN blast_results br ON fs.flanking_id = br.flanking_id
    LEFT JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
    LEFT JOIN annotations a ON br.blast_result_id = a.blast_result_id
    ", where_clause
  )
  
  # Execute query with proper parameter handling
  num_placeholders <- nchar(loci_summary_query) - nchar(gsub("\\?", "", loci_summary_query))
  if (verbose) message("    - Loci summary query placeholders: ", num_placeholders, ", params length: ", length(params))
  
  if (num_placeholders > 0 && length(params) == num_placeholders) {
    loci_summary <- DBI::dbGetQuery(con, loci_summary_query, params)
  } else if (num_placeholders == 0 && length(params) == 0) {
    loci_summary <- DBI::dbGetQuery(con, loci_summary_query)
  } else {
    stop("Parameter mismatch in loci_summary_query: ", num_placeholders, " placeholders, ", length(params), " parameters")
  }
  
  # 2. GO profile (using enhanced function with candidate loci)
  if (verbose) message("  - Analyzing GO term profile...")
  if (verbose) message("    - Calling summarize_go_annotations with candidate_loci rows: ", 
                      if (!is.null(candidate_loci)) nrow(candidate_loci) else "NULL")
  go_profile <- summarize_go_annotations(con, blast_param_id = blast_param_id, 
                                        candidate_loci = candidate_loci, verbose = verbose)
  
  # 3. Protein profile (using existing function)
  if (verbose) message("  - Analyzing protein profile...")
  protein_profile <- summarize_uniprot_annotations(con, blast_param_id = blast_param_id, verbose = FALSE)
  
  # 4. KEGG pathway profile (if requested and available)
  pathway_profile <- NULL
  if (include_pathways && "kegg_references" %in% tables) {
    if (verbose) message("  - Analyzing KEGG pathway profile...")
    
    pathway_query <- paste0("
      SELECT 
        kr.kegg_id,
        kr.pathway_name,
        COUNT(kr.kegg_ref_id) as frequency,
        COUNT(DISTINCT a.annotation_id) as unique_annotations,
        COUNT(DISTINCT br.flanking_id) as unique_loci
      FROM kegg_references kr
      JOIN annotations a ON kr.annotation_id = a.annotation_id
      JOIN blast_results br ON a.blast_result_id = br.blast_result_id
      JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
      JOIN flanking_sequences fs ON br.flanking_id = fs.flanking_id
      JOIN vcf_data vd ON fs.vcf_id = vd.vcf_id
      ", where_clause, "
      GROUP BY kr.kegg_id, kr.pathway_name
      ORDER BY frequency DESC
      LIMIT 50
    ")
    
    # Execute pathway query with proper parameter handling
    num_placeholders <- lengths(regmatches(pathway_query, gregexpr("\\?", pathway_query)))
    if (num_placeholders > 0 && length(params) > 0) {
      pathway_profile <- DBI::dbGetQuery(con, pathway_query, params)
    } else {
      pathway_profile <- DBI::dbGetQuery(con, pathway_query)
    }
  }
  
  # 5. BLAST quality statistics (if requested)
  blast_quality <- NULL
  if (include_blast_stats) {
    if (verbose) message("  - Computing BLAST quality statistics...")
    
    blast_quality_query <- paste0("
      SELECT 
        COUNT(br.blast_result_id) as total_blast_hits,
        COUNT(DISTINCT br.flanking_id) as loci_with_hits,
        ROUND(AVG(br.percent_identity), 2) as avg_percent_identity,
        ROUND(AVG(br.alignment_length), 2) as avg_alignment_length,
        ROUND(AVG(br.e_value), 8) as avg_e_value,
        ROUND(AVG(br.bit_score), 2) as avg_bit_score,
        MIN(br.e_value) as best_e_value,
        MAX(br.bit_score) as best_bit_score,
        COUNT(CASE WHEN br.percent_identity >= 90 THEN 1 END) as high_identity_hits,
        COUNT(CASE WHEN br.e_value <= 1e-10 THEN 1 END) as high_significance_hits
      FROM blast_results br
      JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
      JOIN flanking_sequences fs ON br.flanking_id = fs.flanking_id
      JOIN vcf_data vd ON fs.vcf_id = vd.vcf_id
      ", where_clause
    )
    
    # Execute BLAST quality query with proper parameter handling  
    num_placeholders <- lengths(regmatches(blast_quality_query, gregexpr("\\?", blast_quality_query)))
    if (num_placeholders > 0 && length(params) > 0) {
      blast_quality <- DBI::dbGetQuery(con, blast_quality_query, params)
    } else {
      blast_quality <- DBI::dbGetQuery(con, blast_quality_query)
    }
  }
  
  # 6. Functional diversity metrics
  if (verbose) message("  - Computing functional diversity metrics...")
  
  diversity_query <- paste0("
    WITH functional_stats AS (
      SELECT 
        COUNT(DISTINCT gt.go_id) as unique_go_terms,
        COUNT(DISTINCT gt.go_category) as go_categories,
        COUNT(DISTINCT a.uniprot_accession) as unique_proteins,
        COUNT(DISTINCT CASE WHEN a.gene_names IS NOT NULL AND a.gene_names != '' THEN 
                TRIM(SUBSTR(a.gene_names, 1, CASE 
                  WHEN INSTR(a.gene_names, ';') > 0 THEN INSTR(a.gene_names, ';') - 1 
                  ELSE LENGTH(a.gene_names) 
                END)) END) as unique_gene_names,
        COUNT(DISTINCT br.flanking_id) as annotated_loci
      FROM blast_results br
      JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
      JOIN flanking_sequences fs ON br.flanking_id = fs.flanking_id
      JOIN vcf_data vd ON fs.vcf_id = vd.vcf_id
      LEFT JOIN annotations a ON br.blast_result_id = a.blast_result_id
      LEFT JOIN go_terms gt ON a.annotation_id = gt.annotation_id
      ", where_clause, "
    )
    SELECT 
      unique_go_terms,
      go_categories,
      unique_proteins,
      unique_gene_names,
      annotated_loci,
      ROUND(unique_go_terms * 1.0 / NULLIF(annotated_loci, 0), 2) as go_terms_per_locus,
      ROUND(unique_proteins * 1.0 / NULLIF(annotated_loci, 0), 2) as proteins_per_locus
    FROM functional_stats
  ")
  
  # Execute diversity query with proper parameter handling
  num_placeholders <- lengths(regmatches(diversity_query, gregexpr("\\?", diversity_query)))
  if (num_placeholders > 0 && length(params) > 0) {
    functional_diversity <- DBI::dbGetQuery(con, diversity_query, params)
  } else {
    functional_diversity <- DBI::dbGetQuery(con, diversity_query)
  }
  
  # Clean up temporary table if created
  if (!is.null(candidate_loci)) {
    DBI::dbExecute(con, paste0("DROP TABLE ", temp_table_name))
  }
  
  # Compile results
  result <- list(
    loci_summary = loci_summary,
    go_profile = go_profile,
    protein_profile = protein_profile,
    pathway_profile = pathway_profile,
    blast_quality = blast_quality,
    functional_diversity = functional_diversity,
    parameters = list(
      has_candidate_loci = !is.null(candidate_loci),
      n_candidate_loci = if (!is.null(candidate_loci)) nrow(candidate_loci) else NULL,
      blast_param_id = blast_param_id,
      include_blast_stats = include_blast_stats,
      include_pathways = include_pathways
    )
  )
  
  if (verbose) {
    message("Functional profile complete:")
    if (nrow(loci_summary) > 0) {
      message("  - Input loci: ", loci_summary$total_input_loci[1])
      message("  - Annotated loci: ", loci_summary$annotated_loci[1], 
              " (", loci_summary$annotation_coverage_percent[1], "%)")
      message("  - Total annotations: ", loci_summary$total_annotations[1])
      message("  - Unique proteins: ", loci_summary$unique_proteins[1])
    }
  }
  
  return(result)
}

#' Get annotation coverage statistics across the database
#'
#' Provides detailed statistics about annotation coverage across different
#' annotation types, BLAST parameters, and genomic regions.
#'
#' @param con Database connection object
#' @param by_chromosome Logical. Include chromosome-wise breakdown. Default is TRUE.
#' @param by_blast_param Logical. Include BLAST parameter breakdown. Default is TRUE.
#' @param verbose Logical. Print progress information. Default is TRUE.
#'
#' @return List containing coverage statistics:
#' \itemize{
#'   \item overall: Overall coverage across all annotation types
#'   \item by_type: Coverage broken down by annotation type (GO, UniProt, KEGG)
#'   \item by_chromosome: Coverage by chromosome (if by_chromosome = TRUE)
#'   \item by_blast_param: Coverage by BLAST parameter set (if by_blast_param = TRUE)
#'   \item quality_metrics: Annotation quality indicators
#' }
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#' 
#' # Get comprehensive coverage statistics
#' coverage <- get_annotation_coverage(con)
#' print(coverage$overall)
#' 
#' # Get coverage without chromosome breakdown (faster)
#' coverage_simple <- get_annotation_coverage(con, by_chromosome = FALSE)
#' 
#' close_funseq_db(con)
#' }
#'
#' @export
get_annotation_coverage <- function(con, by_chromosome = TRUE, by_blast_param = TRUE, verbose = TRUE) {
  
  if (verbose) message("Computing annotation coverage statistics...")
  
  # Check if required tables exist
  tables <- DBI::dbListTables(con)
  required_tables <- c("vcf_data", "flanking_sequences", "blast_results", "annotations")
  missing_tables <- required_tables[!required_tables %in% tables]
  
  if (length(missing_tables) > 0) {
    stop("Missing required tables: ", paste(missing_tables, collapse = ", "))
  }
  
  # 1. Overall coverage statistics
  if (verbose) message("  - Computing overall coverage...")
  
  overall_query <- "
    SELECT 
      COUNT(DISTINCT vd.vcf_id) as total_loci,
      COUNT(DISTINCT fs.flanking_id) as loci_with_sequences,
      COUNT(DISTINCT CASE WHEN br.blast_result_id IS NOT NULL THEN fs.flanking_id END) as loci_with_blast_hits,
      COUNT(DISTINCT CASE WHEN a.annotation_id IS NOT NULL THEN fs.flanking_id END) as loci_with_annotations,
      COUNT(DISTINCT CASE WHEN gt.go_term_id IS NOT NULL THEN fs.flanking_id END) as loci_with_go_terms,
      COUNT(DISTINCT CASE WHEN kr.kegg_ref_id IS NOT NULL THEN fs.flanking_id END) as loci_with_kegg,
      ROUND(COUNT(DISTINCT fs.flanking_id) * 100.0 / COUNT(DISTINCT vd.vcf_id), 2) as sequence_coverage_percent,
      ROUND(COUNT(DISTINCT CASE WHEN br.blast_result_id IS NOT NULL THEN fs.flanking_id END) * 100.0 / 
            COUNT(DISTINCT fs.flanking_id), 2) as blast_hit_percent,
      ROUND(COUNT(DISTINCT CASE WHEN a.annotation_id IS NOT NULL THEN fs.flanking_id END) * 100.0 / 
            COUNT(DISTINCT fs.flanking_id), 2) as annotation_percent,
      ROUND(COUNT(DISTINCT CASE WHEN gt.go_term_id IS NOT NULL THEN fs.flanking_id END) * 100.0 / 
            COUNT(DISTINCT fs.flanking_id), 2) as go_coverage_percent,
      ROUND(COUNT(DISTINCT CASE WHEN kr.kegg_ref_id IS NOT NULL THEN fs.flanking_id END) * 100.0 / 
            COUNT(DISTINCT fs.flanking_id), 2) as kegg_coverage_percent
    FROM vcf_data vd
    LEFT JOIN flanking_sequences fs ON vd.vcf_id = fs.vcf_id
    LEFT JOIN blast_results br ON fs.flanking_id = br.flanking_id
    LEFT JOIN annotations a ON br.blast_result_id = a.blast_result_id
    LEFT JOIN go_terms gt ON a.annotation_id = gt.annotation_id
    LEFT JOIN kegg_references kr ON a.annotation_id = kr.annotation_id
  "
  
  overall <- DBI::dbGetQuery(con, overall_query)
  
  # 2. Coverage by annotation type
  if (verbose) message("  - Computing coverage by annotation type...")
  
  by_type_query <- "
    SELECT 
      'GO Terms' as annotation_type,
      COUNT(DISTINCT gt.go_id) as unique_terms,
      COUNT(gt.go_term_id) as total_annotations,
      COUNT(DISTINCT a.annotation_id) as annotated_sequences,
      COUNT(DISTINCT fs.flanking_id) as annotated_loci
    FROM go_terms gt
    JOIN annotations a ON gt.annotation_id = a.annotation_id
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
    JOIN flanking_sequences fs ON br.flanking_id = fs.flanking_id
    
    UNION ALL
    
    SELECT 
      'UniProt Proteins' as annotation_type,
      COUNT(DISTINCT a.uniprot_accession) as unique_terms,
      COUNT(a.annotation_id) as total_annotations,
      COUNT(DISTINCT a.annotation_id) as annotated_sequences,
      COUNT(DISTINCT fs.flanking_id) as annotated_loci
    FROM annotations a
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
    JOIN flanking_sequences fs ON br.flanking_id = fs.flanking_id
    
    UNION ALL
    
    SELECT 
      'KEGG Pathways' as annotation_type,
      COUNT(DISTINCT kr.kegg_id) as unique_terms,
      COUNT(kr.kegg_ref_id) as total_annotations,
      COUNT(DISTINCT a.annotation_id) as annotated_sequences,
      COUNT(DISTINCT fs.flanking_id) as annotated_loci
    FROM kegg_references kr
    JOIN annotations a ON kr.annotation_id = a.annotation_id
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
    JOIN flanking_sequences fs ON br.flanking_id = fs.flanking_id
    
    ORDER BY annotated_loci DESC
  "
  
  by_type <- DBI::dbGetQuery(con, by_type_query)
  
  # 3. Coverage by chromosome (if requested)
  by_chromosome_result <- NULL
  if (by_chromosome) {
    if (verbose) message("  - Computing coverage by chromosome...")
    
    by_chromosome_query <- "
      SELECT 
        vd.chromosome,
        COUNT(DISTINCT vd.vcf_id) as total_loci,
        COUNT(DISTINCT fs.flanking_id) as loci_with_sequences,
        COUNT(DISTINCT CASE WHEN br.blast_result_id IS NOT NULL THEN fs.flanking_id END) as loci_with_blast_hits,
        COUNT(DISTINCT CASE WHEN a.annotation_id IS NOT NULL THEN fs.flanking_id END) as loci_with_annotations,
        ROUND(COUNT(DISTINCT CASE WHEN a.annotation_id IS NOT NULL THEN fs.flanking_id END) * 100.0 / 
              COUNT(DISTINCT fs.flanking_id), 2) as annotation_percent
      FROM vcf_data vd
      LEFT JOIN flanking_sequences fs ON vd.vcf_id = fs.vcf_id
      LEFT JOIN blast_results br ON fs.flanking_id = br.flanking_id
      LEFT JOIN annotations a ON br.blast_result_id = a.blast_result_id
      GROUP BY vd.chromosome
      ORDER BY total_loci DESC
    "
    
    by_chromosome_result <- DBI::dbGetQuery(con, by_chromosome_query)
  }
  
  # 4. Coverage by BLAST parameter (if requested)
  by_blast_param_result <- NULL
  if (by_blast_param) {
    if (verbose) message("  - Computing coverage by BLAST parameter...")
    
    by_blast_param_query <- "
      SELECT 
        bp.blast_param_id,
        bp.blast_type,
        bp.db_name,
        bp.e_value,
        COUNT(DISTINCT br.flanking_id) as loci_with_hits,
        COUNT(DISTINCT a.annotation_id) as total_annotations,
        COUNT(DISTINCT a.uniprot_accession) as unique_proteins,
        COUNT(DISTINCT gt.go_id) as unique_go_terms,
        ROUND(AVG(br.percent_identity), 2) as avg_percent_identity,
        ROUND(AVG(br.e_value), 8) as avg_e_value
      FROM blast_parameters bp
      LEFT JOIN blast_results br ON bp.blast_param_id = br.blast_param_id
      LEFT JOIN annotations a ON br.blast_result_id = a.blast_result_id
      LEFT JOIN go_terms gt ON a.annotation_id = gt.annotation_id
      GROUP BY bp.blast_param_id, bp.blast_type, bp.db_name, bp.e_value
      ORDER BY bp.blast_param_id
    "
    
    by_blast_param_result <- DBI::dbGetQuery(con, by_blast_param_query)
  }
  
  # 5. Quality metrics
  if (verbose) message("  - Computing quality metrics...")
  
  quality_query <- "
    WITH annotation_completeness AS (
      SELECT 
        fs.flanking_id,
        COUNT(DISTINCT a.annotation_id) as num_annotations,
        COUNT(DISTINCT CASE WHEN a.gene_names IS NOT NULL AND a.gene_names != '' THEN a.annotation_id END) as annotations_with_genes,
        COUNT(DISTINCT gt.go_id) as num_go_terms,
        COUNT(DISTINCT gt.go_category) as num_go_categories,
        COUNT(DISTINCT kr.kegg_id) as num_kegg_pathways
      FROM flanking_sequences fs
      LEFT JOIN blast_results br ON fs.flanking_id = br.flanking_id
      LEFT JOIN annotations a ON br.blast_result_id = a.blast_result_id
      LEFT JOIN go_terms gt ON a.annotation_id = gt.annotation_id
      LEFT JOIN kegg_references kr ON a.annotation_id = kr.annotation_id
      GROUP BY fs.flanking_id
    )
    SELECT 
      COUNT(CASE WHEN num_annotations > 0 THEN 1 END) as loci_with_any_annotation,
      COUNT(CASE WHEN num_annotations >= 2 THEN 1 END) as loci_with_multiple_annotations,
      COUNT(CASE WHEN annotations_with_genes > 0 THEN 1 END) as loci_with_gene_names,
      COUNT(CASE WHEN num_go_terms >= 3 THEN 1 END) as loci_with_rich_go_annotation,
      COUNT(CASE WHEN num_go_categories >= 2 THEN 1 END) as loci_with_multiple_go_categories,
      COUNT(CASE WHEN num_kegg_pathways > 0 THEN 1 END) as loci_with_pathway_annotation,
      ROUND(AVG(num_annotations), 2) as avg_annotations_per_locus,
      ROUND(AVG(num_go_terms), 2) as avg_go_terms_per_locus,
      MAX(num_annotations) as max_annotations_per_locus,
      MAX(num_go_terms) as max_go_terms_per_locus
    FROM annotation_completeness
  "
  
  quality_metrics <- DBI::dbGetQuery(con, quality_query)
  
  # Compile results
  result <- list(
    overall = overall,
    by_type = by_type,
    by_chromosome = by_chromosome_result,
    by_blast_param = by_blast_param_result,
    quality_metrics = quality_metrics,
    parameters = list(
      by_chromosome = by_chromosome,
      by_blast_param = by_blast_param
    )
  )
  
  if (verbose) {
    message("Coverage analysis complete:")
    if (nrow(overall) > 0) {
      message("  - Total loci: ", overall$total_loci[1])
      message("  - Annotation coverage: ", overall$annotation_percent[1], "%")
      message("  - GO coverage: ", overall$go_coverage_percent[1], "%")
    }
  }
  
  return(result)
}

#' Set up eggNOG database integration
#'
#' Configures eggNOG database connections and creates necessary tables for
#' storing eggNOG orthologous group annotations, COG categories, and KEGG pathways.
#'
#' @param con Database connection object
#' @param eggnog_db_path Character. Path to local eggNOG database files (optional)
#' @param eggnog_version Character. Version of eggNOG to use. Default is "5.0"
#' @param create_tables Logical. Create eggNOG annotation tables. Default is TRUE
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return List containing setup information:
#' \itemize{
#'   \item tables_created: Names of tables created
#'   \item eggnog_config: Configuration parameters
#'   \item setup_date: Date/time of setup
#' }
#'
#' @details
#' This function sets up the infrastructure for eggNOG annotations by creating
#' the necessary database tables and configuration entries. eggNOG provides
#' orthologous groups, COG functional categories, and KEGG pathway annotations.
#'
#' The following tables are created:
#' \itemize{
#'   \item eggnog_annotations: Links proteins to eggNOG orthologous groups
#'   \item eggnog_ortholog_groups: eggNOG group information and descriptions
#'   \item cog_categories: COG functional category assignments
#'   \item eggnog_kegg_pathways: KEGG pathway annotations from eggNOG
#' }
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#' 
#' # Basic setup with default parameters
#' eggnog_setup <- setup_eggnog_database(con)
#' 
#' # Setup with custom eggNOG database path
#' eggnog_setup <- setup_eggnog_database(
#'   con, 
#'   eggnog_db_path = "/path/to/eggnog/data",
#'   eggnog_version = "5.0"
#' )
#' 
#' close_funseq_db(con)
#' }
#'
#' @export
setup_eggnog_database <- function(con, eggnog_db_path = NULL, eggnog_version = "5.0", 
                                 create_tables = TRUE, verbose = TRUE) {
  
  if (verbose) message("Setting up eggNOG database integration...")
  
  tables_created <- c()
  
  if (create_tables) {
    if (verbose) message("  - Creating eggNOG annotation tables...")
    
    # 1. eggNOG annotations table - links proteins to orthologous groups
    if (verbose) message("    - Creating eggnog_annotations table...")
    DBI::dbExecute(con, "
      CREATE TABLE IF NOT EXISTS eggnog_annotations (
        eggnog_annotation_id INTEGER PRIMARY KEY,
        annotation_id INTEGER NOT NULL,
        eggnog_og_id TEXT NOT NULL,
        og_name TEXT,
        og_description TEXT,
        taxonomic_level TEXT,
        annotation_score REAL,
        annotation_date TEXT NOT NULL,
        FOREIGN KEY (annotation_id) REFERENCES annotations (annotation_id)
      )
    ")
    tables_created <- c(tables_created, "eggnog_annotations")
    
    # 2. eggNOG orthologous groups table - group definitions and metadata
    if (verbose) message("    - Creating eggnog_ortholog_groups table...")
    DBI::dbExecute(con, "
      CREATE TABLE IF NOT EXISTS eggnog_ortholog_groups (
        og_id TEXT PRIMARY KEY,
        og_name TEXT,
        og_description TEXT,
        taxonomic_level TEXT NOT NULL,
        protein_count INTEGER,
        species_count INTEGER,
        eggnog_version TEXT,
        creation_date TEXT
      )
    ")
    tables_created <- c(tables_created, "eggnog_ortholog_groups")
    
    # 3. COG categories table - functional category assignments
    if (verbose) message("    - Creating cog_categories table...")
    DBI::dbExecute(con, "
      CREATE TABLE IF NOT EXISTS cog_categories (
        cog_category_id INTEGER PRIMARY KEY,
        eggnog_annotation_id INTEGER NOT NULL,
        cog_letter TEXT NOT NULL,
        cog_category TEXT NOT NULL,
        cog_description TEXT,
        FOREIGN KEY (eggnog_annotation_id) REFERENCES eggnog_annotations (eggnog_annotation_id)
      )
    ")
    tables_created <- c(tables_created, "cog_categories")
    
    # 4. eggNOG KEGG pathways table - pathway annotations from eggNOG
    if (verbose) message("    - Creating eggnog_kegg_pathways table...")
    DBI::dbExecute(con, "
      CREATE TABLE IF NOT EXISTS eggnog_kegg_pathways (
        eggnog_kegg_id INTEGER PRIMARY KEY,
        eggnog_annotation_id INTEGER NOT NULL,
        kegg_pathway_id TEXT NOT NULL,
        pathway_name TEXT,
        pathway_description TEXT,
        FOREIGN KEY (eggnog_annotation_id) REFERENCES eggnog_annotations (eggnog_annotation_id)
      )
    ")
    tables_created <- c(tables_created, "eggnog_kegg_pathways")
    
    # Create indexes for better performance
    if (verbose) message("    - Creating indexes...")
    DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_eggnog_annotation_id ON eggnog_annotations (annotation_id)")
    DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_eggnog_og_id ON eggnog_annotations (eggnog_og_id)")
    DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_cog_categories_annotation ON cog_categories (eggnog_annotation_id)")
    DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_cog_letter ON cog_categories (cog_letter)")
    DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_eggnog_kegg_annotation ON eggnog_kegg_pathways (eggnog_annotation_id)")
    DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_eggnog_kegg_pathway ON eggnog_kegg_pathways (kegg_pathway_id)")
  }
  
  # Store configuration in metadata table
  if (verbose) message("  - Storing eggNOG configuration...")
  
  config <- list(
    eggnog_version = eggnog_version,
    eggnog_db_path = eggnog_db_path,
    setup_date = Sys.time(),
    tables_created = tables_created
  )
  
  # Store configuration in metadata
  config_json <- jsonlite::toJSON(config, auto_unbox = TRUE)
  
  # Update or insert eggNOG configuration
  DBI::dbExecute(con, "
    INSERT OR REPLACE INTO metadata (key, value) VALUES ('eggnog_config', ?)
  ", list(config_json))
  
  if (verbose) {
    message("eggNOG database setup complete:")
    message("  - Version: ", eggnog_version)
    message("  - Tables created: ", length(tables_created))
    message("  - Database ready for eggNOG annotations")
  }
  
  return(config)
}

#' Annotate proteins with eggNOG orthologous groups
#'
#' Adds eggNOG orthologous group annotations to existing protein annotations,
#' including COG functional categories and KEGG pathway information.
#'
#' @param con Database connection object
#' @param eggnog_results_file Character. Path to eggNOG-mapper results file
#' @param blast_param_id Integer. BLAST parameter set to annotate. If NULL, annotates all.
#' @param overwrite_existing Logical. Overwrite existing eggNOG annotations. Default is FALSE
#' @param min_score Numeric. Minimum annotation score threshold. Default is 0
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return List containing annotation results:
#' \itemize{
#'   \item annotations_added: Number of new eggNOG annotations
#'   \item ortholog_groups_added: Number of new orthologous groups
#'   \item cog_categories_added: Number of COG category assignments
#'   \item kegg_pathways_added: Number of KEGG pathway assignments
#'   \item processing_summary: Summary of processing statistics
#' }
#'
#' @details
#' This function processes eggNOG-mapper output and integrates it with existing
#' protein annotations in the database. eggNOG provides:
#' 
#' \itemize{
#'   \item Orthologous Groups (OGs): Functionally related protein families
#'   \item COG Categories: Functional classification system
#'   \item KEGG Pathways: Metabolic and regulatory pathway annotations
#'   \item Taxonomic Levels: Evolutionary context of annotations
#' }
#'
#' The eggNOG-mapper results file should be in standard tab-delimited format
#' with columns for query protein, eggNOG OG, description, taxonomic level,
#' COG categories, and KEGG pathways.
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#' 
#' # First setup eggNOG database tables
#' setup_eggnog_database(con)
#' 
#' # Annotate with eggNOG results
#' eggnog_results <- annotate_with_eggnog(
#'   con, 
#'   eggnog_results_file = "eggnog_mapper_results.emapper.annotations",
#'   blast_param_id = 1
#' )
#' 
#' print(eggnog_results$processing_summary)
#' 
#' close_funseq_db(con)
#' }
#'
#' @export
annotate_with_eggnog <- function(con, eggnog_results_file, blast_param_id = NULL, 
                                overwrite_existing = FALSE, min_score = 0, verbose = TRUE) {
  
  if (verbose) message("Adding eggNOG annotations from: ", basename(eggnog_results_file))
  
  # Check if file exists
  if (!file.exists(eggnog_results_file)) {
    stop("eggNOG results file not found: ", eggnog_results_file)
  }
  
  # Check if eggNOG tables exist
  tables <- DBI::dbListTables(con)
  required_tables <- c("eggnog_annotations", "eggnog_ortholog_groups", "cog_categories", "eggnog_kegg_pathways")
  missing_tables <- required_tables[!required_tables %in% tables]
  
  if (length(missing_tables) > 0) {
    message("Missing eggNOG tables: ", paste(missing_tables, collapse = ", "))
    message("Running setup_eggnog_database() first...")
    setup_eggnog_database(con, verbose = FALSE)
  }
  
  # Read eggNOG results file
  if (verbose) message("  - Reading eggNOG results file...")
  
  # eggNOG-mapper output format (skip header lines starting with #)
  eggnog_data <- tryCatch({
    # Read all lines and filter out comments
    all_lines <- readLines(eggnog_results_file)
    data_lines <- all_lines[!grepl("^#", all_lines)]
    
    if (length(data_lines) == 0) {
      stop("No data lines found in eggNOG results file")
    }
    
    # Parse tab-delimited data
    # Standard eggNOG-mapper columns: query, seed_ortholog, evalue, score, eggNOG_OGs, max_annot_lvl, COG_category, Description, Preferred_name, GOs, EC, KEGG_ko, KEGG_Pathway, KEGG_Module, KEGG_Reaction, KEGG_rclass, BRITE, KEGG_TC, CAZy, BiGG_Reaction, PFAMs
    
    # Create a connection to read data
    temp_file <- tempfile(fileext = ".tsv")
    writeLines(data_lines, temp_file)
    
    data <- read.table(temp_file, sep = "\t", header = FALSE, quote = "", 
                      stringsAsFactors = FALSE, fill = TRUE, comment.char = "")
    
    # Set column names based on standard eggNOG-mapper output
    expected_cols <- c("query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", 
                      "max_annot_lvl", "COG_category", "Description", "Preferred_name",
                      "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", 
                      "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", 
                      "BiGG_Reaction", "PFAMs")
    
    # Use available columns
    colnames(data) <- expected_cols[1:ncol(data)]
    
    unlink(temp_file)
    data
    
  }, error = function(e) {
    stop("Error reading eggNOG results file: ", e$message)
  })
  
  if (verbose) message("  - Loaded ", nrow(eggnog_data), " eggNOG annotation records")
  
  # Filter by score if specified
  if (min_score > 0) {
    if ("score" %in% colnames(eggnog_data)) {
      eggnog_data <- eggnog_data[!is.na(eggnog_data$score) & eggnog_data$score >= min_score, ]
      if (verbose) message("  - Filtered to ", nrow(eggnog_data), " records with score >= ", min_score)
    }
  }
  
  # Get existing protein annotations to match with eggNOG results
  if (verbose) message("  - Matching with existing protein annotations...")
  
  base_where <- ""
  params <- list()
  
  if (!is.null(blast_param_id)) {
    base_where <- "WHERE bp.blast_param_id = ?"
    params <- list(blast_param_id)
    if (verbose) message("  - Filtering for blast_param_id: ", blast_param_id)
  }
  
  protein_query <- paste0("
    SELECT 
      a.annotation_id,
      a.uniprot_accession,
      a.entry_name,
      br.blast_param_id
    FROM annotations a
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
    JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
    ", base_where, "
    ORDER BY a.annotation_id
  ")
  
  existing_proteins <- DBI::dbGetQuery(con, protein_query, params)
  
  if (verbose) message("  - Found ", nrow(existing_proteins), " existing protein annotations")
  
  # Match eggNOG results with existing annotations
  # Try matching by UniProt accession first, then by entry name
  matched_annotations <- c()
  annotations_added <- 0
  ortholog_groups_added <- 0
  cog_categories_added <- 0
  kegg_pathways_added <- 0
  
  if (verbose) message("  - Processing eggNOG annotations...")
  
  # Track processing statistics
  total_processed <- 0
  matched_count <- 0
  skipped_existing <- 0
  
  for (i in 1:nrow(eggnog_data)) {
    total_processed <- total_processed + 1
    
    if (total_processed %% 1000 == 0 && verbose) {
      message("    - Processed ", total_processed, "/", nrow(eggnog_data), " records...")
    }
    
    query_id <- eggnog_data$query[i]
    
    # Try to match with existing annotations
    # First try exact UniProt accession match
    protein_match <- existing_proteins[existing_proteins$uniprot_accession == query_id, ]
    
    # If no match, try entry name
    if (nrow(protein_match) == 0) {
      protein_match <- existing_proteins[existing_proteins$entry_name == query_id, ]
    }
    
    # If still no match, try partial matches
    if (nrow(protein_match) == 0) {
      # Try matching query ID as substring of UniProt accession
      protein_match <- existing_proteins[grepl(query_id, existing_proteins$uniprot_accession, fixed = TRUE), ]
    }
    
    if (nrow(protein_match) == 0) {
      next  # Skip if no match found
    }
    
    matched_count <- matched_count + 1
    annotation_id <- protein_match$annotation_id[1]  # Use first match if multiple
    
    # Check if eggNOG annotation already exists
    if (!overwrite_existing) {
      existing_check <- DBI::dbGetQuery(con, "
        SELECT COUNT(*) as count FROM eggnog_annotations WHERE annotation_id = ?
      ", list(annotation_id))
      
      if (existing_check$count[1] > 0) {
        skipped_existing <- skipped_existing + 1
        next
      }
    }
    
    # Extract eggNOG annotation data
    og_id <- eggnog_data$eggNOG_OGs[i]
    description <- if ("Description" %in% colnames(eggnog_data)) eggnog_data$Description[i] else ""
    preferred_name <- if ("Preferred_name" %in% colnames(eggnog_data)) eggnog_data$Preferred_name[i] else ""
    max_annot_lvl <- if ("max_annot_lvl" %in% colnames(eggnog_data)) eggnog_data$max_annot_lvl[i] else ""
    score <- if ("score" %in% colnames(eggnog_data)) eggnog_data$score[i] else NA
    
    # Skip if no eggNOG OG assignment
    if (is.na(og_id) || og_id == "" || og_id == "-") {
      next
    }
    
    # Add eggNOG annotation
    DBI::dbExecute(con, "
      INSERT OR REPLACE INTO eggnog_annotations 
      (annotation_id, eggnog_og_id, og_name, og_description, taxonomic_level, annotation_score, annotation_date)
      VALUES (?, ?, ?, ?, ?, ?, ?)
    ", list(annotation_id, og_id, preferred_name, description, max_annot_lvl, score, Sys.time()))
    
    annotations_added <- annotations_added + 1
    eggnog_annotation_id <- DBI::dbGetQuery(con, "SELECT last_insert_rowid() as id")$id[1]
    
    # Add orthologous group information
    DBI::dbExecute(con, "
      INSERT OR IGNORE INTO eggnog_ortholog_groups 
      (og_id, og_name, og_description, taxonomic_level, eggnog_version, creation_date)
      VALUES (?, ?, ?, ?, ?, ?)
    ", list(og_id, preferred_name, description, max_annot_lvl, "5.0", Sys.time()))
    
    if (DBI::dbGetQuery(con, "SELECT changes() as changes")$changes[1] > 0) {
      ortholog_groups_added <- ortholog_groups_added + 1
    }
    
    # Add COG categories
    if ("COG_category" %in% colnames(eggnog_data)) {
      cog_categories <- eggnog_data$COG_category[i]
      if (!is.na(cog_categories) && cog_categories != "" && cog_categories != "-") {
        # COG categories are often single letters or combinations
        cog_letters <- strsplit(cog_categories, "")[[1]]
        
        # COG category descriptions (simplified mapping)
        cog_descriptions <- c(
          "A" = "RNA processing and modification",
          "B" = "Chromatin structure and dynamics", 
          "C" = "Energy production and conversion",
          "D" = "Cell cycle control, cell division, chromosome partitioning",
          "E" = "Amino acid transport and metabolism",
          "F" = "Nucleotide transport and metabolism",
          "G" = "Carbohydrate transport and metabolism",
          "H" = "Coenzyme transport and metabolism",
          "I" = "Lipid transport and metabolism",
          "J" = "Translation, ribosomal structure and biogenesis",
          "K" = "Transcription",
          "L" = "Replication, recombination and repair",
          "M" = "Cell wall/membrane/envelope biogenesis",
          "N" = "Cell motility",
          "O" = "Posttranslational modification, protein turnover, chaperones",
          "P" = "Inorganic ion transport and metabolism",
          "Q" = "Secondary metabolites biosynthesis, transport and catabolism",
          "R" = "General function prediction only",
          "S" = "Function unknown",
          "T" = "Signal transduction mechanisms",
          "U" = "Intracellular trafficking, secretion, and vesicular transport",
          "V" = "Defense mechanisms",
          "W" = "Extracellular structures",
          "X" = "Mobilome: prophages, transposons",
          "Y" = "Nuclear structure",
          "Z" = "Cytoskeleton"
        )
        
        for (cog_letter in cog_letters) {
          if (cog_letter %in% names(cog_descriptions)) {
            DBI::dbExecute(con, "
              INSERT INTO cog_categories 
              (eggnog_annotation_id, cog_letter, cog_category, cog_description)
              VALUES (?, ?, ?, ?)
            ", list(eggnog_annotation_id, cog_letter, cog_letter, cog_descriptions[cog_letter]))
            
            cog_categories_added <- cog_categories_added + 1
          }
        }
      }
    }
    
    # Add KEGG pathways
    if ("KEGG_Pathway" %in% colnames(eggnog_data)) {
      kegg_pathways <- eggnog_data$KEGG_Pathway[i]
      if (!is.na(kegg_pathways) && kegg_pathways != "" && kegg_pathways != "-") {
        # KEGG pathways are often comma-separated
        pathway_list <- strsplit(kegg_pathways, ",")[[1]]
        
        for (pathway in pathway_list) {
          pathway <- trimws(pathway)
          if (pathway != "") {
            # Extract pathway ID and name if formatted as "ko12345:pathway_name"
            if (grepl(":", pathway)) {
              parts <- strsplit(pathway, ":")[[1]]
              pathway_id <- parts[1]
              pathway_name <- if (length(parts) > 1) parts[2] else ""
            } else {
              pathway_id <- pathway
              pathway_name <- ""
            }
            
            DBI::dbExecute(con, "
              INSERT INTO eggnog_kegg_pathways 
              (eggnog_annotation_id, kegg_pathway_id, pathway_name, pathway_description)
              VALUES (?, ?, ?, ?)
            ", list(eggnog_annotation_id, pathway_id, pathway_name, ""))
            
            kegg_pathways_added <- kegg_pathways_added + 1
          }
        }
      }
    }
  }
  
  # Compile results
  processing_summary <- data.frame(
    total_eggnog_records = nrow(eggnog_data),
    total_processed = total_processed,
    matched_with_proteins = matched_count,
    skipped_existing = skipped_existing,
    annotations_added = annotations_added,
    ortholog_groups_added = ortholog_groups_added,
    cog_categories_added = cog_categories_added,
    kegg_pathways_added = kegg_pathways_added,
    processing_date = Sys.time()
  )
  
  result <- list(
    annotations_added = annotations_added,
    ortholog_groups_added = ortholog_groups_added,
    cog_categories_added = cog_categories_added,
    kegg_pathways_added = kegg_pathways_added,
    processing_summary = processing_summary
  )
  
  if (verbose) {
    message("eggNOG annotation complete:")
    message("  - Annotations added: ", annotations_added)
    message("  - Ortholog groups: ", ortholog_groups_added)
    message("  - COG categories: ", cog_categories_added)
    message("  - KEGG pathways: ", kegg_pathways_added)
    message("  - Match rate: ", round(matched_count/total_processed*100, 1), "%")
  }
  
  return(result)
}

#' Summarize COG categories from eggNOG annotations
#'
#' Provides comprehensive analysis of COG (Clusters of Orthologous Groups) functional 
#' categories from eggNOG annotations, including distribution patterns and coverage statistics.
#'
#' @param con Database connection object
#' @param blast_param_id Integer. Specific BLAST parameter set to analyze. If NULL, analyzes all.
#' @param include_functional_groups Logical. Group COG categories by major functional areas. Default is TRUE.
#' @param min_frequency Integer. Minimum frequency threshold for COG categories to include. Default is 1.
#' @param verbose Logical. Print progress information. Default is TRUE.
#'
#' @return List containing COG category analysis:
#' \itemize{
#'   \item overview: Summary statistics of COG category coverage
#'   \item category_distribution: Frequency of each COG category
#'   \item functional_groups: Analysis by major functional areas (if include_functional_groups = TRUE)
#'   \item annotation_coverage: Coverage statistics across annotations
#'   \item category_descriptions: Full descriptions of COG categories found
#' }
#'
#' @details
#' COG categories are organized into major functional groups:
#' \itemize{
#'   \item Information Storage and Processing: Categories A, B, J, K, L
#'   \item Cellular Processes and Signaling: Categories D, M, N, O, T, U, V, W, Y, Z
#'   \item Metabolism: Categories C, E, F, G, H, I, P, Q
#'   \item Poorly Characterized: Categories R, S
#' }
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#' 
#' # Analyze all COG categories
#' cog_summary <- summarize_cog_categories(con)
#' print(cog_summary$overview)
#' print(cog_summary$functional_groups)
#' 
#' # Focus on specific BLAST parameter set
#' cog_summary_blast1 <- summarize_cog_categories(con, blast_param_id = 1)
#' 
#' close_funseq_db(con)
#' }
#'
#' @export
summarize_cog_categories <- function(con, blast_param_id = NULL, include_functional_groups = TRUE, 
                                   min_frequency = 1, verbose = TRUE) {
  
  if (verbose) message("Summarizing COG categories from eggNOG annotations...")
  
  # Check if required tables exist
  tables <- DBI::dbListTables(con)
  required_tables <- c("eggnog_annotations", "cog_categories", "annotations", "blast_results")
  missing_tables <- required_tables[!required_tables %in% tables]
  
  if (length(missing_tables) > 0) {
    stop("Missing required tables for COG analysis: ", paste(missing_tables, collapse = ", "),
         "\nPlease run setup_eggnog_database() and annotate_with_eggnog() first.")
  }
  
  # Build base query with optional blast_param_id filter
  base_where <- ""
  params <- list()
  
  if (!is.null(blast_param_id)) {
    base_where <- "WHERE bp.blast_param_id = ?"
    params <- list(blast_param_id)
    if (verbose) message("  - Filtering for blast_param_id: ", blast_param_id)
  }
  
  # 1. Overview statistics
  if (verbose) message("  - Computing COG category overview...")
  
  overview_query <- paste0("
    SELECT 
      COUNT(DISTINCT cc.cog_letter) as unique_cog_categories,
      COUNT(cc.cog_category_id) as total_cog_assignments,
      COUNT(DISTINCT ea.eggnog_annotation_id) as proteins_with_cog,
      COUNT(DISTINCT a.annotation_id) as annotated_sequences,
      COUNT(DISTINCT br.flanking_id) as annotated_loci
    FROM cog_categories cc
    JOIN eggnog_annotations ea ON cc.eggnog_annotation_id = ea.eggnog_annotation_id
    JOIN annotations a ON ea.annotation_id = a.annotation_id
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
    JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
    ", base_where
  )
  
  if (length(params) > 0) {
    overview <- DBI::dbGetQuery(con, overview_query, params)
  } else {
    overview <- DBI::dbGetQuery(con, overview_query)
  }
  
  # 2. COG category distribution
  if (verbose) message("  - Analyzing COG category distribution...")
  
  distribution_query <- paste0("
    SELECT 
      cc.cog_letter,
      cc.cog_category,
      cc.cog_description,
      COUNT(cc.cog_category_id) as frequency,
      COUNT(DISTINCT ea.eggnog_annotation_id) as unique_proteins,
      COUNT(DISTINCT br.flanking_id) as unique_loci
    FROM cog_categories cc
    JOIN eggnog_annotations ea ON cc.eggnog_annotation_id = ea.eggnog_annotation_id
    JOIN annotations a ON ea.annotation_id = a.annotation_id
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
    JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
    ", base_where, "
    GROUP BY cc.cog_letter, cc.cog_category, cc.cog_description
    HAVING COUNT(cc.cog_category_id) >= ?
    ORDER BY frequency DESC
  ")
  
  distribution_params <- c(params, list(min_frequency))
  # Execute distribution query with proper parameter handling
  num_placeholders <- nchar(distribution_query) - nchar(gsub("\\?", "", distribution_query))
  if (num_placeholders > 0 && length(distribution_params) == num_placeholders) {
    category_distribution <- DBI::dbGetQuery(con, distribution_query, distribution_params)
  } else if (num_placeholders == 0 && length(distribution_params) == 0) {
    category_distribution <- DBI::dbGetQuery(con, distribution_query)
  } else {
    stop("Parameter mismatch in distribution_query: ", num_placeholders, " placeholders, ", length(distribution_params), " parameters")
  }
  
  # 3. Functional group analysis (if requested)
  functional_groups <- NULL
  if (include_functional_groups) {
    if (verbose) message("  - Analyzing functional group patterns...")
    
    # Define COG functional groups
    functional_group_mapping <- data.frame(
      cog_letter = c("A", "B", "J", "K", "L",  # Information Storage and Processing
                     "D", "M", "N", "O", "T", "U", "V", "W", "Y", "Z",  # Cellular Processes
                     "C", "E", "F", "G", "H", "I", "P", "Q",  # Metabolism
                     "R", "S"),  # Poorly Characterized
      functional_group = c(rep("Information Storage and Processing", 5),
                          rep("Cellular Processes and Signaling", 10),
                          rep("Metabolism", 8),
                          rep("Poorly Characterized", 2)),
      stringsAsFactors = FALSE
    )
    
    if (nrow(category_distribution) > 0) {
      # Add functional group information to category distribution
      distribution_with_groups <- merge(category_distribution, functional_group_mapping, 
                                      by = "cog_letter", all.x = TRUE)
      
      # Summarize by functional group
      functional_groups <- aggregate(
        cbind(frequency, unique_proteins, unique_loci) ~ functional_group,
        data = distribution_with_groups,
        FUN = sum,
        na.rm = TRUE
      )
      
      # Add percentage calculations
      total_assignments <- sum(functional_groups$frequency)
      functional_groups$percent_of_assignments <- round(
        (functional_groups$frequency / total_assignments) * 100, 2
      )
      
      # Order by frequency
      functional_groups <- functional_groups[order(functional_groups$frequency, decreasing = TRUE), ]
    } else {
      functional_groups <- data.frame(
        functional_group = character(0),
        frequency = integer(0),
        unique_proteins = integer(0),
        unique_loci = integer(0),
        percent_of_assignments = numeric(0)
      )
    }
  }
  
  # 4. Annotation coverage statistics
  if (verbose) message("  - Computing annotation coverage...")
  
  coverage_query <- paste0("
    WITH cog_coverage AS (
      SELECT 
        bp.blast_param_id,
        bp.blast_type,
        bp.db_name,
        COUNT(DISTINCT cc.cog_letter) as unique_cog_categories,
        COUNT(cc.cog_category_id) as total_cog_assignments,
        COUNT(DISTINCT ea.eggnog_annotation_id) as proteins_with_cog,
        COUNT(DISTINCT br.flanking_id) as loci_with_cog
      FROM blast_parameters bp
      LEFT JOIN blast_results br ON bp.blast_param_id = br.blast_param_id
      LEFT JOIN annotations a ON br.blast_result_id = a.blast_result_id
      LEFT JOIN eggnog_annotations ea ON a.annotation_id = ea.annotation_id
      LEFT JOIN cog_categories cc ON ea.eggnog_annotation_id = cc.eggnog_annotation_id
      ", if (base_where == "") "" else gsub("WHERE", "AND", base_where), "
      GROUP BY bp.blast_param_id, bp.blast_type, bp.db_name
      ORDER BY bp.blast_param_id
    )
    SELECT * FROM cog_coverage
  ")
  
  if (length(params) > 0 && !is.null(blast_param_id)) {
    # For specific blast_param_id, add comparison with overall
    coverage_params <- params
  } else {
    coverage_params <- list()
  }
  
  # Execute coverage query with proper parameter handling
  num_placeholders <- nchar(coverage_query) - nchar(gsub("\\?", "", coverage_query))
  if (num_placeholders > 0 && length(coverage_params) == num_placeholders) {
    annotation_coverage <- DBI::dbGetQuery(con, coverage_query, coverage_params)
  } else if (num_placeholders == 0 && length(coverage_params) == 0) {
    annotation_coverage <- DBI::dbGetQuery(con, coverage_query)
  } else {
    stop("Parameter mismatch in coverage_query: ", num_placeholders, " placeholders, ", length(coverage_params), " parameters")
  }
  
  # 5. Create comprehensive category descriptions
  category_descriptions <- if (nrow(category_distribution) > 0) {
    unique(category_distribution[, c("cog_letter", "cog_category", "cog_description")])
  } else {
    data.frame(
      cog_letter = character(0),
      cog_category = character(0), 
      cog_description = character(0)
    )
  }
  
  # Compile results
  result <- list(
    overview = overview,
    category_distribution = category_distribution,
    functional_groups = functional_groups,
    annotation_coverage = annotation_coverage,
    category_descriptions = category_descriptions,
    parameters = list(
      blast_param_id = blast_param_id,
      include_functional_groups = include_functional_groups,
      min_frequency = min_frequency
    )
  )
  
  if (verbose) {
    message("COG category analysis complete:")
    if (nrow(overview) > 0 && overview$total_cog_assignments[1] > 0) {
      message("  - Unique COG categories: ", overview$unique_cog_categories[1])
      message("  - Total COG assignments: ", overview$total_cog_assignments[1])
      message("  - Proteins with COG: ", overview$proteins_with_cog[1])
      message("  - Annotated loci: ", overview$annotated_loci[1])
      
      if (include_functional_groups && !is.null(functional_groups) && nrow(functional_groups) > 0) {
        message("  - Top functional group: ", 
                functional_groups$functional_group[1], " (", 
                functional_groups$percent_of_assignments[1], "% of assignments)")
      }
    } else {
      message("  - No COG annotations found with current filters")
    }
  }
  
  return(result)
}

# INTERNAL HELPER FUNCTIONS

#' Create GO functional modules from GO terms
#'
#' Groups GO terms into biologically meaningful functional modules based on 
#' semantic relationships and common functional themes.
#'
#' @param go_terms_data Data frame with go_id, go_term, and go_category columns
#' @param verbose Logical. Print progress information. Default is TRUE.
#'
#' @return List containing functional modules with their associated GO terms
#'
.create_go_functional_modules <- function(go_terms_data, verbose = TRUE) {
  
  if (is.null(go_terms_data) || nrow(go_terms_data) == 0) {
    return(list())
  }
  
  # Define functional module mappings based on GO term patterns
  module_patterns <- list(
    "DNA Repair and Maintenance" = c(
      "DNA repair", "DNA damage", "DNA recombination", "DNA replication", 
      "chromosome organization", "chromatin", "nucleotide excision repair",
      "base excision repair", "mismatch repair", "double-strand break repair"
    ),
    
    "Cell Cycle and Division" = c(
      "cell cycle", "cell division", "mitotic", "meiotic", "cytokinesis",
      "chromosome segregation", "spindle", "centrosome", "kinetochore"
    ),
    
    "Transcription and Gene Expression" = c(
      "transcription", "RNA polymerase", "promoter", "enhancer", "silencer",
      "chromatin remodeling", "histone", "gene expression", "RNA processing"
    ),
    
    "Translation and Protein Synthesis" = c(
      "translation", "ribosome", "tRNA", "aminoacyl", "protein synthesis",
      "ribosomal", "peptide bond", "codon"
    ),
    
    "Protein Folding and Quality Control" = c(
      "protein folding", "chaperone", "heat shock", "unfolded protein",
      "endoplasmic reticulum stress", "proteasome", "ubiquitin", "protein degradation"
    ),
    
    "Metabolic Processes" = c(
      "metabolic process", "glycolysis", "gluconeogenesis", "citrate cycle",
      "oxidative phosphorylation", "fatty acid", "amino acid metabolism",
      "carbohydrate metabolism", "lipid metabolism"
    ),
    
    "Signal Transduction" = c(
      "signal transduction", "receptor", "kinase", "phosphatase", "second messenger",
      "G protein", "calcium signaling", "cAMP", "phosphorylation"
    ),
    
    "Transport and Membrane" = c(
      "transport", "membrane", "ion channel", "carrier", "ATPase",
      "vesicle", "endocytosis", "exocytosis", "secretion"
    ),
    
    "Immune Response" = c(
      "immune", "antigen", "antibody", "T cell", "B cell", "cytokine",
      "complement", "inflammation", "defense response"
    ),
    
    "Development and Differentiation" = c(
      "development", "differentiation", "morphogenesis", "pattern specification",
      "embryonic", "cell fate", "organogenesis", "tissue development"
    ),
    
    "Apoptosis and Cell Death" = c(
      "apoptosis", "cell death", "programmed cell death", "caspase",
      "cytochrome c", "Bcl-2", "p53", "DNA fragmentation"
    ),
    
    "Oxidative Stress Response" = c(
      "oxidative stress", "antioxidant", "reactive oxygen", "superoxide",
      "catalase", "peroxidase", "glutathione", "oxidoreductase"
    )
  )
  
  # Initialize results structure
  modules <- list()
  
  for (module_name in names(module_patterns)) {
    patterns <- module_patterns[[module_name]]
    
    # Find GO terms matching any of the patterns for this module
    matching_terms <- go_terms_data[
      apply(go_terms_data, 1, function(row) {
        term_text <- tolower(paste(row["go_term"], collapse = " "))
        any(sapply(patterns, function(pattern) grepl(pattern, term_text, ignore.case = TRUE)))
      }), 
    ]
    
    if (nrow(matching_terms) > 0) {
      modules[[module_name]] <- list(
        terms = matching_terms,
        term_count = nrow(matching_terms),
        categories = unique(matching_terms$go_category),
        total_frequency = sum(matching_terms$frequency, na.rm = TRUE)
      )
    }
  }
  
  # Add unclassified terms
  classified_terms <- unique(do.call(rbind, lapply(modules, function(m) m$terms))$go_id)
  unclassified <- go_terms_data[!go_terms_data$go_id %in% classified_terms, ]
  
  if (nrow(unclassified) > 0) {
    modules[["Other/Unclassified"]] <- list(
      terms = unclassified,
      term_count = nrow(unclassified),
      categories = unique(unclassified$go_category),
      total_frequency = sum(unclassified$frequency, na.rm = TRUE)
    )
  }
  
  # Sort modules by total frequency
  modules <- modules[order(sapply(modules, function(m) m$total_frequency), decreasing = TRUE)]
  
  if (verbose) {
    message("    - Created ", length(modules), " functional modules")
    message("    - Classified ", sum(sapply(modules, function(m) m$term_count)), " GO terms")
  }
  
  return(modules)
}

#' Create KEGG functional modules from pathway data
#'
#' Groups KEGG pathways into biologically meaningful functional modules based on 
#' KEGG pathway classification and functional categories.
#'
#' @param pathway_data Data frame with kegg_id, pathway_name, and frequency columns
#' @param verbose Logical. Print progress information. Default is TRUE.
#'
#' @return List containing functional modules with their associated KEGG pathways
#'
.create_kegg_functional_modules <- function(pathway_data, verbose = TRUE) {
  
  if (is.null(pathway_data) || nrow(pathway_data) == 0) {
    return(list())
  }
  
  # Define KEGG functional module mappings based on pathway patterns
  module_patterns <- list(
    "Carbohydrate Metabolism" = c(
      "glycolysis", "gluconeogenesis", "pentose phosphate", "citrate cycle",
      "pyruvate metabolism", "glyoxylate", "propanoate", "butanoate",
      "starch", "sucrose", "galactose", "fructose", "mannose"
    ),
    
    "Energy Metabolism" = c(
      "oxidative phosphorylation", "photosynthesis", "carbon fixation",
      "methane metabolism", "nitrogen metabolism", "sulfur metabolism"
    ),
    
    "Lipid Metabolism" = c(
      "fatty acid", "steroid", "glycerolipid", "glycerophospholipid",
      "sphingolipid", "arachidonic acid", "linoleic acid", "bile acid"
    ),
    
    "Amino Acid Metabolism" = c(
      "alanine", "arginine", "aspartate", "cysteine", "glutamate", 
      "glycine", "histidine", "isoleucine", "leucine", "lysine",
      "methionine", "phenylalanine", "proline", "serine", "threonine",
      "tryptophan", "tyrosine", "valine"
    ),
    
    "Nucleotide Metabolism" = c(
      "purine metabolism", "pyrimidine metabolism", "nucleotide sugar"
    ),
    
    "Vitamin and Cofactor Metabolism" = c(
      "thiamine", "riboflavin", "niacin", "pantothenate", "pyridoxine",
      "biotin", "folate", "ubiquinone", "terpenoid", "porphyrin"
    ),
    
    "Signal Transduction" = c(
      "MAPK", "calcium", "phosphatidylinositol", "phospholipase",
      "sphingolipid signaling", "mTOR", "AMPK", "insulin", "Wnt", "Notch", "Hedgehog"
    ),
    
    "Cell Cycle and Proliferation" = c(
      "cell cycle", "DNA replication", "mismatch repair", "homologous recombination",
      "non-homologous end-joining", "p53", "cellular senescence"
    ),
    
    "Membrane Transport" = c(
      "ABC transporters", "phosphotransferase", "bacterial secretion",
      "protein export", "sulfur relay"
    ),
    
    "Immune System" = c(
      "complement", "toll-like receptor", "NOD-like receptor", "cytosolic DNA-sensing",
      "natural killer", "T cell receptor", "B cell receptor", "chemokine",
      "IL-17", "Th1", "Th2", "Th17", "Treg"
    ),
    
    "Endocrine System" = c(
      "insulin", "adipocytokine", "PPAR", "thyroid hormone", "estrogen",
      "androgen", "progesterone", "cortisol", "aldosterone", "parathyroid"
    ),
    
    "Development and Regeneration" = c(
      "axon guidance", "osteoclast", "melanogenesis", "adipogenesis",
      "chondroitin sulfate", "heparan sulfate"
    ),
    
    "Drug Metabolism" = c(
      "cytochrome P450", "drug metabolism", "xenobiotics"
    ),
    
    "Disease Pathways" = c(
      "cancer", "alzheimer", "parkinson", "huntington", "diabetes",
      "pathways in cancer", "viral carcinogenesis", "chemical carcinogenesis"
    ),
    
    # Organism-specific gene modules (for cases where pathway names aren't available)
    "Zebrafish Genes" = c(
      "danio rerio", "zebrafish", "dre:"
    ),
    
    "Fugu Genes" = c(
      "takifugu rubripes", "fugu", "tru:"
    ),
    
    "Human Genes" = c(
      "homo sapiens", "human", "hsa:"
    ),
    
    "Mouse Genes" = c(
      "mus musculus", "mouse", "mmu:"
    ),
    
    "Rat Genes" = c(
      "rattus norvegicus", "rat", "rno:"
    )
  )
  
  # Initialize results structure
  modules <- list()
  
  for (module_name in names(module_patterns)) {
    patterns <- module_patterns[[module_name]]
    
    # Find pathways matching any of the patterns for this module
    matching_pathways <- pathway_data[
      apply(pathway_data, 1, function(row) {
        # Check both kegg_id and pathway_name for matches
        pathway_text <- tolower(paste(row["kegg_id"], row["pathway_name"], collapse = " "))
        any(sapply(patterns, function(pattern) grepl(pattern, pathway_text, ignore.case = TRUE)))
      }), 
    ]
    
    if (nrow(matching_pathways) > 0) {
      modules[[module_name]] <- list(
        pathways = matching_pathways,
        pathway_count = nrow(matching_pathways),
        total_frequency = sum(matching_pathways$frequency, na.rm = TRUE),
        unique_loci = sum(matching_pathways$unique_loci, na.rm = TRUE)
      )
    }
  }
  
  # Add unclassified pathways
  classified_pathways <- unique(do.call(rbind, lapply(modules, function(m) m$pathways))$kegg_id)
  unclassified <- pathway_data[!pathway_data$kegg_id %in% classified_pathways, ]
  
  if (nrow(unclassified) > 0) {
    modules[["Other/Unclassified"]] <- list(
      pathways = unclassified,
      pathway_count = nrow(unclassified),
      total_frequency = sum(unclassified$frequency, na.rm = TRUE),
      unique_loci = sum(unclassified$unique_loci, na.rm = TRUE)
    )
  }
  
  # Sort modules by total frequency
  modules <- modules[order(sapply(modules, function(m) m$total_frequency), decreasing = TRUE)]
  
  if (verbose) {
    message("    - Created ", length(modules), " pathway modules")
    message("    - Classified ", sum(sapply(modules, function(m) m$pathway_count)), " pathways")
  }
  
  return(modules)
}

#' Create integrated functional summary from multiple annotation types
#'
#' Integrates functional modules from GO, KEGG, and COG to identify common themes
#' and create a unified functional characterization.
#'
#' @param go_modules List. GO functional modules from .create_go_functional_modules()
#' @param kegg_modules List. KEGG pathway modules from .create_kegg_functional_modules()
#' @param cog_modules Data frame. COG functional groups 
#' @param verbose Logical. Print progress information. Default is TRUE.
#'
#' @return List containing integrated functional themes
#'
.create_integrated_functional_summary <- function(go_modules = NULL, kegg_modules = NULL, 
                                                 cog_modules = NULL, verbose = TRUE) {
  
  # Define functional theme mappings that integrate across annotation types
  integrated_themes <- list()
  
  # Theme 1: DNA and Genome Maintenance
  dna_theme <- list(
    name = "DNA and Genome Maintenance",
    description = "Processes involved in DNA repair, replication, and genome stability",
    sources = list()
  )
  
  if (!is.null(go_modules)) {
    dna_go <- go_modules[grepl("DNA Repair|Cell Cycle", names(go_modules))]
    if (length(dna_go) > 0) {
      dna_theme$sources$GO <- dna_go
    }
  }
  
  if (!is.null(kegg_modules)) {
    dna_kegg <- kegg_modules[grepl("Cell Cycle|DNA", names(kegg_modules))]
    if (length(dna_kegg) > 0) {
      dna_theme$sources$KEGG <- dna_kegg
    }
  }
  
  if (!is.null(cog_modules)) {
    dna_cog <- cog_modules[grepl("Replication|Recombination", cog_modules$functional_group), ]
    if (nrow(dna_cog) > 0) {
      dna_theme$sources$COG <- dna_cog
    }
  }
  
  if (length(dna_theme$sources) > 0) {
    integrated_themes[["DNA and Genome Maintenance"]] <- dna_theme
  }
  
  # Theme 2: Metabolism and Energy
  metabolism_theme <- list(
    name = "Metabolism and Energy",
    description = "Metabolic processes and energy production pathways",
    sources = list()
  )
  
  if (!is.null(go_modules)) {
    metab_go <- go_modules[grepl("Metabolic", names(go_modules))]
    if (length(metab_go) > 0) {
      metabolism_theme$sources$GO <- metab_go
    }
  }
  
  if (!is.null(kegg_modules)) {
    metab_kegg <- kegg_modules[grepl("Metabolism|Energy", names(kegg_modules))]
    if (length(metab_kegg) > 0) {
      metabolism_theme$sources$KEGG <- metab_kegg
    }
  }
  
  if (!is.null(cog_modules)) {
    metab_cog <- cog_modules[grepl("Energy|Amino acid|Carbohydrate|Lipid", cog_modules$functional_group), ]
    if (nrow(metab_cog) > 0) {
      metabolism_theme$sources$COG <- metab_cog
    }
  }
  
  if (length(metabolism_theme$sources) > 0) {
    integrated_themes[["Metabolism and Energy"]] <- metabolism_theme
  }
  
  # Theme 3: Protein Processing and Quality Control
  protein_theme <- list(
    name = "Protein Processing and Quality Control",
    description = "Protein synthesis, folding, modification, and degradation",
    sources = list()
  )
  
  if (!is.null(go_modules)) {
    prot_go <- go_modules[grepl("Translation|Protein Folding", names(go_modules))]
    if (length(prot_go) > 0) {
      protein_theme$sources$GO <- prot_go
    }
  }
  
  if (!is.null(cog_modules)) {
    prot_cog <- cog_modules[grepl("Translation|Posttranslational", cog_modules$functional_group), ]
    if (nrow(prot_cog) > 0) {
      protein_theme$sources$COG <- prot_cog
    }
  }
  
  if (length(protein_theme$sources) > 0) {
    integrated_themes[["Protein Processing and Quality Control"]] <- protein_theme
  }
  
  # Theme 4: Signaling and Regulation
  signaling_theme <- list(
    name = "Signaling and Regulation",
    description = "Signal transduction and regulatory mechanisms",
    sources = list()
  )
  
  if (!is.null(go_modules)) {
    sig_go <- go_modules[grepl("Signal Transduction|Transcription", names(go_modules))]
    if (length(sig_go) > 0) {
      signaling_theme$sources$GO <- sig_go
    }
  }
  
  if (!is.null(kegg_modules)) {
    sig_kegg <- kegg_modules[grepl("Signal Transduction", names(kegg_modules))]
    if (length(sig_kegg) > 0) {
      signaling_theme$sources$KEGG <- sig_kegg
    }
  }
  
  if (!is.null(cog_modules)) {
    sig_cog <- cog_modules[grepl("Signal transduction|Transcription", cog_modules$functional_group), ]
    if (nrow(sig_cog) > 0) {
      signaling_theme$sources$COG <- sig_cog
    }
  }
  
  if (length(signaling_theme$sources) > 0) {
    integrated_themes[["Signaling and Regulation"]] <- signaling_theme
  }
  
  # Theme 5: Stress Response and Defense
  stress_theme <- list(
    name = "Stress Response and Defense",
    description = "Cellular stress response and defense mechanisms",
    sources = list()
  )
  
  if (!is.null(go_modules)) {
    stress_go <- go_modules[grepl("Oxidative Stress|Immune|Apoptosis", names(go_modules))]
    if (length(stress_go) > 0) {
      stress_theme$sources$GO <- stress_go
    }
  }
  
  if (!is.null(kegg_modules)) {
    stress_kegg <- kegg_modules[grepl("Immune", names(kegg_modules))]
    if (length(stress_kegg) > 0) {
      stress_theme$sources$KEGG <- stress_kegg
    }
  }
  
  if (!is.null(cog_modules)) {
    stress_cog <- cog_modules[grepl("Defense", cog_modules$functional_group), ]
    if (nrow(stress_cog) > 0) {
      stress_theme$sources$COG <- stress_cog
    }
  }
  
  if (length(stress_theme$sources) > 0) {
    integrated_themes[["Stress Response and Defense"]] <- stress_theme
  }
  
  if (verbose && length(integrated_themes) > 0) {
    message("    - Created ", length(integrated_themes), " integrated functional themes")
  }
  
  return(integrated_themes)
}

#' Analyze overlap between functional modules from different annotation types
#'
#' Identifies common functional themes and overlapping coverage between
#' GO, KEGG, and COG annotation systems.
#'
#' @param go_modules List. GO functional modules
#' @param kegg_modules List. KEGG pathway modules
#' @param cog_modules Data frame. COG functional groups
#' @param verbose Logical. Print progress information. Default is TRUE.
#'
#' @return List containing overlap analysis results
#'
.analyze_module_overlap <- function(go_modules = NULL, kegg_modules = NULL, 
                                   cog_modules = NULL, verbose = TRUE) {
  
  overlap_results <- list()
  
  # Count modules per annotation type
  module_counts <- list(
    GO = if (!is.null(go_modules)) length(go_modules) else 0,
    KEGG = if (!is.null(kegg_modules)) length(kegg_modules) else 0,
    COG = if (!is.null(cog_modules)) nrow(cog_modules) else 0
  )
  
  overlap_results$module_counts <- module_counts
  
  # Identify common functional themes by name similarity
  all_module_names <- c()
  
  if (!is.null(go_modules)) {
    all_module_names <- c(all_module_names, paste("GO:", names(go_modules)))
  }
  
  if (!is.null(kegg_modules)) {
    all_module_names <- c(all_module_names, paste("KEGG:", names(kegg_modules)))
  }
  
  if (!is.null(cog_modules)) {
    all_module_names <- c(all_module_names, paste("COG:", cog_modules$functional_group))
  }
  
  # Find potential overlaps based on keyword matching
  common_keywords <- c("metabolism", "transport", "signal", "protein", "DNA", "cell cycle", "immune")
  
  keyword_overlap <- list()
  for (keyword in common_keywords) {
    matching_modules <- all_module_names[grepl(keyword, all_module_names, ignore.case = TRUE)]
    if (length(matching_modules) > 1) {
      keyword_overlap[[keyword]] <- matching_modules
    }
  }
  
  overlap_results$keyword_overlap <- keyword_overlap
  
  # Calculate coverage statistics
  coverage_stats <- list()
  
  if (!is.null(go_modules) && !is.null(kegg_modules)) {
    # Count how many themes are covered by both GO and KEGG
    go_themes <- tolower(names(go_modules))
    kegg_themes <- tolower(names(kegg_modules))
    
    shared_themes <- sum(sapply(go_themes, function(g) {
      any(sapply(kegg_themes, function(k) {
        length(intersect(strsplit(g, " ")[[1]], strsplit(k, " ")[[1]])) > 0
      }))
    }))
    
    coverage_stats$go_kegg_overlap <- shared_themes
  }
  
  overlap_results$coverage_stats <- coverage_stats
  
  # Summary metrics
  overlap_results$summary <- list(
    total_unique_themes = length(unique(all_module_names)),
    annotation_systems_used = sum(module_counts > 0),
    common_keyword_themes = length(keyword_overlap)
  )
  
  if (verbose) {
    message("    - Analyzed overlap across ", overlap_results$summary$annotation_systems_used, " annotation systems")
    message("    - Found ", overlap_results$summary$common_keyword_themes, " common functional themes")
  }
  
  return(overlap_results)
}

#' Update KEGG pathway names using KEGGREST
#'
#' Fixes missing or incorrect KEGG pathway names in the database by using KEGGREST
#' to fetch actual pathway information from the KEGG API. This addresses the issue
#' where UniProt returns "-" as pathway names.
#'
#' @param con Database connection object
#' @param limit Integer. Maximum number of KEGG IDs to process. Default is NULL (process all)
#' @param batch_size Integer. Number of KEGG IDs to process in each batch. Default is 10
#' @param verbose Logical. Print progress information. Default is TRUE
#' @param dry_run Logical. If TRUE, show what would be updated without making changes. Default is FALSE
#'
#' @return List containing update results:
#' \itemize{
#'   \item total_kegg_ids: Total number of KEGG IDs found in database
#'   \item processed_ids: Number of IDs successfully processed
#'   \item updated_pathways: Number of pathway names updated
#'   \item failed_ids: Vector of KEGG IDs that failed to process
#'   \item processing_summary: Detailed summary of updates
#' }
#'
#' @details
#' This function retrieves KEGG IDs from the database that have missing or invalid
#' pathway names (NULL, empty, or "-"), then uses KEGGREST to fetch the actual
#' pathway information from KEGG. For gene IDs like "dre:563201", it extracts
#' pathway names and updates the database.
#'
#' The function handles different types of KEGG IDs:
#' \itemize{
#'   \item Gene IDs (e.g., "dre:563201"): Fetches gene information and pathway associations
#'   \item Pathway IDs (e.g., "dre00562"): Fetches pathway names directly
#'   \item Module IDs (e.g., "M00130"): Fetches module information
#' }
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#' 
#' # Preview what would be updated (dry run)
#' preview <- update_kegg_pathway_names_with_keggrest(con, limit = 5, dry_run = TRUE)
#' print(preview$processing_summary)
#' 
#' # Actually update pathway names
#' results <- update_kegg_pathway_names_with_keggrest(con, limit = 50)
#' print(results$processing_summary)
#' 
#' close_funseq_db(con)
#' }
#'
#' @export
update_kegg_pathway_names_with_keggrest <- function(con, limit = NULL, batch_size = 10, 
                                                   verbose = TRUE, dry_run = FALSE) {
  
  if (verbose) {
    action <- if (dry_run) "Previewing" else "Updating"
    message(action, " KEGG pathway names using KEGGREST...")
  }
  
  # Check if KEGGREST is available
  if (!requireNamespace("KEGGREST", quietly = TRUE)) {
    stop("KEGGREST package is required but not installed. Please install it with:\n",
         "BiocManager::install('KEGGREST')")
  }
  
  # Get KEGG IDs with missing or invalid pathway names
  if (verbose) message("  - Retrieving KEGG IDs with missing pathway names...")
  
  base_query <- "
    SELECT 
      kr.kegg_ref_id,
      kr.kegg_id,
      kr.pathway_name,
      a.uniprot_accession,
      a.entry_name
    FROM kegg_references kr
    JOIN annotations a ON kr.annotation_id = a.annotation_id
    WHERE kr.pathway_name IS NULL 
       OR kr.pathway_name = '' 
       OR kr.pathway_name = '-'
    ORDER BY kr.kegg_ref_id
  "
  
  if (!is.null(limit)) {
    base_query <- paste(base_query, "LIMIT", limit)
  }
  
  kegg_records <- DBI::dbGetQuery(con, base_query)
  
  if (nrow(kegg_records) == 0) {
    if (verbose) message("  - No KEGG references need updating")
    return(list(
      total_kegg_ids = 0,
      processed_ids = 0,
      updated_pathways = 0,
      failed_ids = character(0),
      processing_summary = data.frame()
    ))
  }
  
  if (verbose) message("  - Found ", nrow(kegg_records), " KEGG references to process")
  
  # Process in batches to avoid overwhelming the API
  total_records <- nrow(kegg_records)
  processed_count <- 0
  updated_count <- 0
  failed_ids <- character(0)
  processing_details <- list()
  
  # Split into batches
  batch_indices <- split(1:total_records, ceiling(seq_along(1:total_records) / batch_size))
  
  for (batch_num in seq_along(batch_indices)) {
    batch_indices_current <- batch_indices[[batch_num]]
    batch_records <- kegg_records[batch_indices_current, ]
    
    if (verbose) {
      message("  - Processing batch ", batch_num, "/", length(batch_indices), 
              " (", nrow(batch_records), " records)...")
    }
    
    for (i in 1:nrow(batch_records)) {
      record <- batch_records[i, ]
      kegg_id <- record$kegg_id
      kegg_ref_id <- record$kegg_ref_id
      
      if (verbose && (processed_count + 1) %% 25 == 0) {
        message("    - Processed ", processed_count + 1, "/", total_records, " records...")
      }
      
      tryCatch({
        # Use KEGGREST to get pathway information
        kegg_data <- KEGGREST::keggGet(kegg_id)
        
        if (length(kegg_data) > 0) {
          gene_info <- kegg_data[[1]]
          new_pathway_names <- character(0)
          
          # Extract pathway information
          if (!is.null(gene_info$PATHWAY)) {
            pathways <- gene_info$PATHWAY
            # Convert to descriptive format: "pathway_id: pathway_name"
            pathway_descriptions <- paste(names(pathways), pathways, sep = ": ")
            new_pathway_names <- pathway_descriptions
          }
          
          # If no pathways, use gene name/symbol information
          if (length(new_pathway_names) == 0) {
            if (!is.null(gene_info$NAME)) {
              new_pathway_names <- paste("Gene:", gene_info$NAME)
            } else if (!is.null(gene_info$SYMBOL)) {
              new_pathway_names <- paste("Gene:", gene_info$SYMBOL)
            } else {
              # Use organism-based naming as fallback
              if (grepl("^[a-z]{3}:", kegg_id)) {
                organism_code <- substr(kegg_id, 1, 3)
                organism_map <- c(
                  "dre" = "Danio rerio (zebrafish)",
                  "tru" = "Takifugu rubripes (fugu)",
                  "hsa" = "Homo sapiens (human)",
                  "mmu" = "Mus musculus (mouse)",
                  "rno" = "Rattus norvegicus (rat)",
                  "ath" = "Arabidopsis thaliana",
                  "sce" = "Saccharomyces cerevisiae (yeast)",
                  "eco" = "Escherichia coli"
                )
                
                if (organism_code %in% names(organism_map)) {
                  new_pathway_names <- paste("Gene from", organism_map[organism_code])
                } else {
                  new_pathway_names <- paste("Gene from organism:", organism_code)
                }
              } else {
                new_pathway_names <- "Unknown pathway"
              }
            }
          }
          
          # Use the first pathway name (or combined if multiple)
          final_pathway_name <- if (length(new_pathway_names) > 1) {
            # If multiple pathways, combine them
            paste(new_pathway_names, collapse = "; ")
          } else {
            new_pathway_names[1]
          }
          
          # Store processing details
          processing_details[[as.character(kegg_ref_id)]] <- list(
            kegg_id = kegg_id,
            old_pathway_name = record$pathway_name,
            new_pathway_name = final_pathway_name,
            pathways_found = length(gene_info$PATHWAY %||% 0),
            status = "success"
          )
          
          # Update database (unless dry run)
          if (!dry_run) {
            DBI::dbExecute(con, "
              UPDATE kegg_references 
              SET pathway_name = ?
              WHERE kegg_ref_id = ?
            ", list(final_pathway_name, kegg_ref_id))
          }
          
          updated_count <- updated_count + 1
          
        } else {
          # No data returned from KEGG
          processing_details[[as.character(kegg_ref_id)]] <- list(
            kegg_id = kegg_id,
            old_pathway_name = record$pathway_name,
            new_pathway_name = NA,
            pathways_found = 0,
            status = "no_data"
          )
        }
        
        processed_count <- processed_count + 1
        
      }, error = function(e) {
        # Handle API errors
        failed_ids <- c(failed_ids, kegg_id)
        processing_details[[as.character(kegg_ref_id)]] <- list(
          kegg_id = kegg_id,
          old_pathway_name = record$pathway_name,
          new_pathway_name = NA,
          pathways_found = 0,
          status = paste("error:", e$message)
        )
        
        if (verbose) {
          message("    - Failed to process ", kegg_id, ": ", e$message)
        }
        
        processed_count <- processed_count + 1
      })
      
      # Rate limiting: small delay between requests
      Sys.sleep(0.1)
    }
    
    # Longer delay between batches
    if (batch_num < length(batch_indices)) {
      Sys.sleep(1)
    }
  }
  
  # Create processing summary
  processing_summary <- do.call(rbind, lapply(names(processing_details), function(ref_id) {
    details <- processing_details[[ref_id]]
    data.frame(
      kegg_ref_id = as.integer(ref_id),
      kegg_id = details$kegg_id,
      old_pathway_name = details$old_pathway_name %||% "",
      new_pathway_name = details$new_pathway_name %||% "",
      pathways_found = details$pathways_found,
      status = details$status,
      stringsAsFactors = FALSE
    )
  }))
  
  # Compile results
  result <- list(
    total_kegg_ids = total_records,
    processed_ids = processed_count,
    updated_pathways = updated_count,
    failed_ids = failed_ids,
    processing_summary = processing_summary
  )
  
  if (verbose) {
    action_past <- if (dry_run) "previewed" else "updated"
    message("KEGG pathway name processing complete:")
    message("  - Total records: ", total_records)
    message("  - Successfully processed: ", processed_count)
    message("  - Pathway names ", action_past, ": ", updated_count)
    if (length(failed_ids) > 0) {
      message("  - Failed IDs: ", length(failed_ids))
    }
    
    if (!dry_run && updated_count > 0) {
      message("  - Database updated with KEGGREST pathway information")
    }
  }
  
  return(result)
}

#' Analyze KEGG modules and pathway networks using ggkegg
#'
#' Performs advanced KEGG pathway analysis using the ggkegg package to extract
#' pathway network information, module relationships, and functional groupings.
#' This provides deeper insight into pathway organization than simple pathway lists.
#'
#' @param con Database connection object
#' @param candidate_loci Data frame with chromosome and position columns, or NULL for all loci
#' @param blast_param_id Integer. Specific BLAST parameter set to analyze. If NULL, analyzes all.
#' @param use_cached_pathways Logical. Use cached ggkegg pathway data. Default is TRUE
#' @param include_network_analysis Logical. Perform pathway network analysis. Default is TRUE
#' @param max_pathways_to_analyze Integer. Maximum number of pathways to analyze with ggkegg. Default is 20
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return List containing KEGG module analysis:
#' \itemize{
#'   \item pathway_summary: Basic pathway frequency statistics
#'   \item module_analysis: KEGG module information and relationships
#'   \item network_data: Pathway network graphs and connectivity (if include_network_analysis = TRUE)
#'   \item functional_modules: Pathways grouped by functional themes
#'   \item pathway_interactions: Cross-pathway gene/compound relationships
#'   \item coverage_stats: Analysis coverage and quality metrics
#' }
#'
#' @details
#' This function extends basic KEGG pathway analysis by using ggkegg to:
#' \itemize{
#'   \item Retrieve pathway network structures from KEGG
#'   \item Identify KEGG modules that group related pathways
#'   \item Analyze gene/compound connectivity between pathways
#'   \item Provide functional module classifications
#'   \item Generate network-based pathway summaries
#' }
#'
#' The analysis leverages both the KEGG database structure and the ggkegg package's
#' ability to parse pathway networks and module relationships.
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#' 
#' # Analyze all KEGG data with module analysis
#' kegg_modules <- analyze_kegg_modules(con)
#' print(kegg_modules$module_analysis)
#' 
#' # Focus on candidate loci with network analysis
#' candidates <- data.frame(chromosome = c("LG1", "LG2"), position = c(12345, 67890))
#' kegg_analysis <- analyze_kegg_modules(
#'   con, 
#'   candidate_loci = candidates,
#'   include_network_analysis = TRUE
#' )
#' 
#' close_funseq_db(con)
#' }
#'
#' @export
analyze_kegg_modules <- function(con, candidate_loci = NULL, blast_param_id = NULL,
                                use_cached_pathways = TRUE, include_network_analysis = TRUE,
                                max_pathways_to_analyze = 20, verbose = TRUE) {
  
  if (verbose) message("Analyzing KEGG modules and pathway networks...")
  
  # Check if ggkegg is available
  ggkegg_available <- requireNamespace("ggkegg", quietly = TRUE)
  if (!ggkegg_available) {
    if (verbose) message("  - ggkegg not available, using basic pathway analysis")
    include_network_analysis <- FALSE
  }
  
  # 1. Get basic pathway data from database
  if (verbose) message("  - Retrieving KEGG pathway data from database...")
  
  # Build base query
  base_where_clause <- "WHERE kr.kegg_id IS NOT NULL"
  params <- list()
  
  # Add candidate loci filter if provided
  if (!is.null(candidate_loci)) {
    if (verbose) message("    - Filtering for candidate loci...")
    if (!all(c("chromosome", "position") %in% colnames(candidate_loci))) {
      stop("candidate_loci must have 'chromosome' and 'position' columns")
    }
    
    # Create temporary table for candidate loci
    temp_table_name <- paste0("temp_candidates_", as.integer(Sys.time()))
    
    # Create and populate temporary table
    DBI::dbExecute(con, paste0("CREATE TEMPORARY TABLE ", temp_table_name, " (chromosome TEXT, position INTEGER)"))
    
    for (i in 1:nrow(candidate_loci)) {
      DBI::dbExecute(con, paste0("INSERT INTO ", temp_table_name, " VALUES (?, ?)"),
                    list(candidate_loci$chromosome[i], candidate_loci$position[i]))
    }
    
    base_where_clause <- paste0(base_where_clause, " AND vd.chromosome IN (SELECT chromosome FROM ", temp_table_name, 
                               ") AND vd.position IN (SELECT position FROM ", temp_table_name, ")")
  }
  
  # Add blast parameter filter if provided
  if (!is.null(blast_param_id)) {
    base_where_clause <- paste0(base_where_clause, " AND bp.blast_param_id = ?")
    params <- c(params, list(blast_param_id))
    if (verbose) message("    - Filtering for blast_param_id: ", blast_param_id)
  }
  
  # Query for pathway data
  pathway_query <- paste0("
    SELECT 
      kr.kegg_id,
      kr.pathway_name,
      COUNT(*) as frequency,
      COUNT(DISTINCT vd.vcf_id) as unique_loci,
      COUNT(DISTINCT a.uniprot_accession) as unique_proteins,
      GROUP_CONCAT(DISTINCT a.uniprot_accession) as protein_list
    FROM kegg_references kr
    JOIN annotations a ON kr.annotation_id = a.annotation_id
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
    JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
    JOIN flanking_sequences fs ON br.flanking_id = fs.flanking_id
    JOIN vcf_data vd ON fs.vcf_id = vd.vcf_id
    ", base_where_clause, "
    GROUP BY kr.kegg_id, kr.pathway_name
    ORDER BY frequency DESC
  ")
  
  # Execute query with proper parameter handling
  num_placeholders <- nchar(pathway_query) - nchar(gsub("\\?", "", pathway_query))
  if (num_placeholders > 0 && length(params) == num_placeholders) {
    pathway_data <- DBI::dbGetQuery(con, pathway_query, params)
  } else if (num_placeholders == 0 && length(params) == 0) {
    pathway_data <- DBI::dbGetQuery(con, pathway_query)
  } else {
    stop("Parameter mismatch in pathway_query: ", num_placeholders, " placeholders, ", length(params), " parameters")
  }
  
  if (nrow(pathway_data) == 0) {
    if (verbose) message("  - No KEGG pathway data found with current filters")
    return(list(
      pathway_summary = data.frame(),
      module_analysis = list(),
      network_data = list(),
      functional_modules = list(),
      pathway_interactions = data.frame(),
      coverage_stats = list()
    ))
  }
  
  if (verbose) message("  - Found ", nrow(pathway_data), " unique KEGG pathways")
  
  # 2. Basic pathway summary
  pathway_summary <- pathway_data[, c("kegg_id", "pathway_name", "frequency", "unique_loci", "unique_proteins")]
  
  # 3. KEGG module analysis using ggkegg (if available)
  module_analysis <- list()
  network_data <- list()
  
  if (ggkegg_available && include_network_analysis) {
    if (verbose) message("  - Performing ggkegg-based module analysis...")
    
    # Focus on most frequent pathways to avoid API overload
    top_pathways <- head(pathway_data, max_pathways_to_analyze)
    
    # Extract pathway IDs that look like pathway maps (e.g., "dre00562")
    pathway_ids <- top_pathways$kegg_id[grepl("^[a-z]{3}\\d{5}$", top_pathways$kegg_id)]
    
    if (length(pathway_ids) > 0) {
      if (verbose) message("    - Analyzing ", length(pathway_ids), " pathway networks...")
      
      pathway_networks <- list()
      module_info <- list()
      
      for (i in seq_along(pathway_ids)) {
        pathway_id <- pathway_ids[i]
        
        if (verbose && i %% 5 == 0) {
          message("      - Processing pathway ", i, "/", length(pathway_ids), ": ", pathway_id)
        }
        
        tryCatch({
          # Get pathway network using ggkegg
          pathway_graph <- ggkegg::pathway(pathway_id, use_cache = use_cached_pathways)
          
          if (!is.null(pathway_graph)) {
            # Extract network information
            nodes <- igraph::V(pathway_graph)
            edges <- igraph::E(pathway_graph)
            
            # Store network data
            pathway_networks[[pathway_id]] <- list(
              pathway_id = pathway_id,
              node_count = length(nodes),
              edge_count = length(edges),
              graph = pathway_graph,
              genes = nodes$name[nodes$type == "gene"],
              compounds = nodes$name[nodes$type == "compound"]
            )
            
            # Try to get associated module information
            # Look for module IDs in the pathway data
            if (!is.null(nodes$pathway_id) && length(unique(nodes$pathway_id)) > 0) {
              # This pathway is part of larger network - could indicate module structure
              module_info[[pathway_id]] <- list(
                pathway_id = pathway_id,
                connected_pathways = unique(nodes$pathway_id),
                module_genes = length(unique(nodes$name[nodes$type == "gene"])),
                module_compounds = length(unique(nodes$name[nodes$type == "compound"]))
              )
            }
          }
          
          # Rate limiting
          Sys.sleep(0.5)
          
        }, error = function(e) {
          if (verbose) {
            message("      - Failed to process pathway ", pathway_id, ": ", e$message)
          }
        })
      }
      
      network_data <- pathway_networks
      module_analysis <- module_info
      
    } else {
      if (verbose) message("    - No valid pathway IDs found for network analysis")
    }
  }
  
  # 4. Functional module classification (using enhanced patterns)
  if (verbose) message("  - Creating functional module classifications...")
  
  functional_modules <- .create_kegg_functional_modules(pathway_data, verbose = verbose)
  
  # 5. Pathway interaction analysis
  pathway_interactions <- data.frame()
  
  if (length(network_data) > 1) {
    if (verbose) message("  - Analyzing pathway interactions...")
    
    # Find shared genes/compounds between pathways
    interaction_data <- list()
    
    pathway_names <- names(network_data)
    for (i in 1:(length(pathway_names) - 1)) {
      for (j in (i + 1):length(pathway_names)) {
        pathway1 <- pathway_names[i]
        pathway2 <- pathway_names[j]
        
        net1 <- network_data[[pathway1]]
        net2 <- network_data[[pathway2]]
        
        shared_genes <- intersect(net1$genes, net2$genes)
        shared_compounds <- intersect(net1$compounds, net2$compounds)
        
        if (length(shared_genes) > 0 || length(shared_compounds) > 0) {
          interaction_data[[paste(pathway1, pathway2, sep = "_")]] <- data.frame(
            pathway1 = pathway1,
            pathway2 = pathway2,
            shared_genes = length(shared_genes),
            shared_compounds = length(shared_compounds),
            total_shared = length(shared_genes) + length(shared_compounds),
            shared_gene_list = paste(shared_genes, collapse = "; "),
            shared_compound_list = paste(shared_compounds, collapse = "; "),
            stringsAsFactors = FALSE
          )
        }
      }
    }
    
    if (length(interaction_data) > 0) {
      pathway_interactions <- do.call(rbind, interaction_data)
      rownames(pathway_interactions) <- NULL
      
      # Sort by total shared elements
      pathway_interactions <- pathway_interactions[order(pathway_interactions$total_shared, decreasing = TRUE), ]
    }
  }
  
  # 6. Coverage statistics
  coverage_stats <- list(
    total_pathways = nrow(pathway_data),
    pathways_with_networks = length(network_data),
    total_modules = length(module_analysis),
    functional_module_count = length(functional_modules),
    pathway_interactions_count = nrow(pathway_interactions),
    analysis_parameters = list(
      ggkegg_available = ggkegg_available,
      network_analysis_performed = include_network_analysis && ggkegg_available,
      max_pathways_analyzed = max_pathways_to_analyze,
      candidate_loci_provided = !is.null(candidate_loci),
      blast_param_filter = blast_param_id
    )
  )
  
  # Compile results
  result <- list(
    pathway_summary = pathway_summary,
    module_analysis = module_analysis,
    network_data = network_data,
    functional_modules = functional_modules,
    pathway_interactions = pathway_interactions,
    coverage_stats = coverage_stats
  )
  
  if (verbose) {
    message("KEGG module analysis complete:")
    message("  - Pathways analyzed: ", coverage_stats$total_pathways)
    if (ggkegg_available && include_network_analysis) {
      message("  - Network data retrieved: ", coverage_stats$pathways_with_networks)
      message("  - Functional modules: ", coverage_stats$functional_module_count)
      message("  - Pathway interactions: ", coverage_stats$pathway_interactions_count)
    } else {
      message("  - Basic functional classification performed")
    }
  }
  
  return(result)
}

#' Create pathway network summary for candidate loci
#'
#' Generates a comprehensive network-based summary of how candidate loci
#' connect through metabolic and regulatory pathways, identifying key
#' pathway hubs and functional relationships.
#'
#' @param con Database connection object
#' @param candidate_loci Data frame with chromosome and position columns
#' @param background_loci Data frame with chromosome and position columns for background, or NULL
#' @param include_pathway_crosstalk Logical. Analyze pathway crosstalk and interactions. Default is TRUE
#' @param include_centrality_analysis Logical. Calculate network centrality metrics. Default is TRUE
#' @param min_pathway_size Integer. Minimum pathway size to include in network. Default is 3
#' @param max_pathways Integer. Maximum pathways to include in network analysis. Default is 50
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return List containing pathway network analysis:
#' \itemize{
#'   \item network_graph: igraph object representing pathway network
#'   \item pathway_hubs: Pathways with high connectivity (centrality)
#'   \item crosstalk_analysis: Pathway interaction analysis
#'   \item candidate_pathway_profile: Pathway enrichment in candidate vs background
#'   \item network_metrics: Network topology metrics
#'   \item functional_clusters: Groups of functionally related pathways
#' }
#'
#' @details
#' This function creates a comprehensive pathway network by:
#' \itemize{
#'   \item Building a network where nodes are pathways and edges represent shared genes/compounds
#'   \item Identifying pathway hubs (highly connected pathways)
#'   \item Analyzing pathway crosstalk and functional relationships
#'   \item Comparing candidate loci pathway profiles to background
#'   \item Clustering pathways into functional modules
#'   \item Computing network topology metrics
#' }
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#' 
#' # Define candidate and background loci
#' candidates <- data.frame(
#'   chromosome = c("LG1", "LG2", "LG3"), 
#'   position = c(12345, 67890, 111213)
#' )
#' background <- data.frame(
#'   chromosome = c("LG1", "LG2", "LG3", "LG4"), 
#'   position = c(54321, 98765, 131415, 161718)
#' )
#' 
#' # Create pathway network summary
#' network_summary <- create_pathway_network_summary(
#'   con, 
#'   candidate_loci = candidates,
#'   background_loci = background
#' )
#' 
#' print(network_summary$pathway_hubs)
#' print(network_summary$network_metrics)
#' 
#' close_funseq_db(con)
#' }
#'
#' @export
create_pathway_network_summary <- function(con, candidate_loci, background_loci = NULL,
                                          include_pathway_crosstalk = TRUE, 
                                          include_centrality_analysis = TRUE,
                                          min_pathway_size = 3, max_pathways = 50,
                                          verbose = TRUE) {
  
  if (verbose) message("Creating pathway network summary for candidate loci...")
  
  # Check if required packages are available
  igraph_available <- requireNamespace("igraph", quietly = TRUE)
  if (!igraph_available) {
    stop("igraph package is required for network analysis. Please install it with:\n",
         "install.packages('igraph')")
  }
  
  # Validate input
  if (is.null(candidate_loci) || nrow(candidate_loci) == 0) {
    stop("candidate_loci is required and must contain at least one locus")
  }
  
  if (!all(c("chromosome", "position") %in% colnames(candidate_loci))) {
    stop("candidate_loci must have 'chromosome' and 'position' columns")
  }
  
  # 1. Get pathway data for candidate loci
  if (verbose) message("  - Retrieving pathway data for candidate loci...")
  
  candidate_pathways <- analyze_kegg_modules(
    con, 
    candidate_loci = candidate_loci,
    include_network_analysis = TRUE,
    max_pathways_to_analyze = max_pathways,
    verbose = FALSE
  )
  
  if (length(candidate_pathways$network_data) == 0) {
    if (verbose) message("  - No pathway network data found for candidate loci")
    return(list(
      network_graph = NULL,
      pathway_hubs = data.frame(),
      crosstalk_analysis = data.frame(),
      candidate_pathway_profile = data.frame(),
      network_metrics = list(),
      functional_clusters = list()
    ))
  }
  
  # 2. Get background pathway data (if provided)
  background_pathways <- NULL
  if (!is.null(background_loci)) {
    if (verbose) message("  - Retrieving pathway data for background loci...")
    
    if (!all(c("chromosome", "position") %in% colnames(background_loci))) {
      stop("background_loci must have 'chromosome' and 'position' columns")
    }
    
    background_pathways <- analyze_kegg_modules(
      con, 
      candidate_loci = background_loci,
      include_network_analysis = FALSE,  # Just need pathway summary
      verbose = FALSE
    )
  }
  
  # 3. Build pathway network
  if (verbose) message("  - Building pathway interaction network...")
  
  # Extract pathway interactions from candidate data
  pathway_interactions <- candidate_pathways$pathway_interactions
  
  if (nrow(pathway_interactions) == 0) {
    if (verbose) message("    - No pathway interactions found")
    network_graph <- NULL
    pathway_hubs <- data.frame()
    network_metrics <- list()
  } else {
    # Create igraph network
    edges <- pathway_interactions[, c("pathway1", "pathway2", "total_shared")]
    colnames(edges) <- c("from", "to", "weight")
    
    # Get all unique pathways
    all_pathways <- unique(c(edges$from, edges$to))
    
    # Filter by minimum pathway size if network data is available
    if (min_pathway_size > 1) {
      pathway_summary <- candidate_pathways$pathway_summary
      valid_pathways <- pathway_summary$kegg_id[pathway_summary$unique_loci >= min_pathway_size]
      
      # Filter edges to only include valid pathways
      edges <- edges[edges$from %in% valid_pathways & edges$to %in% valid_pathways, ]
      all_pathways <- unique(c(edges$from, edges$to))
    }
    
    if (nrow(edges) > 0) {
      # Create network graph
      network_graph <- igraph::graph_from_data_frame(edges, vertices = all_pathways, directed = FALSE)
      
      # Add vertex attributes (pathway information)
      pathway_info <- candidate_pathways$pathway_summary
      vertex_names <- igraph::V(network_graph)$name
      
      for (i in seq_along(vertex_names)) {
        pathway_id <- vertex_names[i]
        pathway_row <- pathway_info[pathway_info$kegg_id == pathway_id, ]
        
        if (nrow(pathway_row) > 0) {
          igraph::V(network_graph)$pathway_name[i] <- pathway_row$pathway_name[1]
          igraph::V(network_graph)$frequency[i] <- pathway_row$frequency[1]
          igraph::V(network_graph)$unique_loci[i] <- pathway_row$unique_loci[1]
        } else {
          igraph::V(network_graph)$pathway_name[i] <- "Unknown"
          igraph::V(network_graph)$frequency[i] <- 0
          igraph::V(network_graph)$unique_loci[i] <- 0
        }
      }
      
      # 4. Calculate centrality measures (if requested)
      pathway_hubs <- data.frame()
      if (include_centrality_analysis) {
        if (verbose) message("  - Calculating pathway centrality measures...")
        
        # Calculate various centrality measures
        degree_centrality <- igraph::degree(network_graph)
        betweenness_centrality <- igraph::betweenness(network_graph)
        closeness_centrality <- igraph::closeness(network_graph)
        eigenvector_centrality <- igraph::eigen_centrality(network_graph)$vector
        
        # Create pathway hubs dataframe
        pathway_hubs <- data.frame(
          pathway_id = names(degree_centrality),
          pathway_name = igraph::V(network_graph)$pathway_name,
          degree = degree_centrality,
          betweenness = betweenness_centrality,
          closeness = closeness_centrality,
          eigenvector = eigenvector_centrality,
          frequency = igraph::V(network_graph)$frequency,
          unique_loci = igraph::V(network_graph)$unique_loci,
          stringsAsFactors = FALSE
        )
        
        # Calculate hub score (composite of centrality measures)
        pathway_hubs$hub_score <- scale(pathway_hubs$degree)[,1] + 
                                 scale(pathway_hubs$betweenness)[,1] + 
                                 scale(pathway_hubs$eigenvector)[,1]
        
        # Sort by hub score
        pathway_hubs <- pathway_hubs[order(pathway_hubs$hub_score, decreasing = TRUE), ]
        rownames(pathway_hubs) <- NULL
      }
      
      # 5. Network topology metrics
      if (verbose) message("  - Computing network topology metrics...")
      
      network_metrics <- list(
        nodes = igraph::vcount(network_graph),
        edges = igraph::ecount(network_graph),
        density = igraph::edge_density(network_graph),
        diameter = ifelse(igraph::is.connected(network_graph), igraph::diameter(network_graph), NA),
        avg_path_length = ifelse(igraph::is.connected(network_graph), igraph::mean_distance(network_graph), NA),
        clustering_coefficient = igraph::transitivity(network_graph),
        connected_components = igraph::components(network_graph)$no,
        largest_component_size = max(igraph::components(network_graph)$csize)
      )
      
    } else {
      network_graph <- NULL
      pathway_hubs <- data.frame()
      network_metrics <- list()
    }
  }
  
  # 6. Pathway crosstalk analysis
  crosstalk_analysis <- data.frame()
  if (include_pathway_crosstalk && !is.null(network_graph)) {
    if (verbose) message("  - Analyzing pathway crosstalk...")
    
    # Use the existing pathway interactions data
    crosstalk_analysis <- pathway_interactions
    
    # Add crosstalk strength categories
    if (nrow(crosstalk_analysis) > 0) {
      crosstalk_analysis$crosstalk_strength <- cut(
        crosstalk_analysis$total_shared,
        breaks = c(0, 1, 3, 10, Inf),
        labels = c("Weak", "Moderate", "Strong", "Very Strong"),
        include.lowest = TRUE
      )
      
      # Sort by total shared elements
      crosstalk_analysis <- crosstalk_analysis[order(crosstalk_analysis$total_shared, decreasing = TRUE), ]
    }
  }
  
  # 7. Candidate vs background pathway profile (if background provided)
  candidate_pathway_profile <- data.frame()
  if (!is.null(background_pathways)) {
    if (verbose) message("  - Comparing candidate vs background pathway profiles...")
    
    candidate_summary <- candidate_pathways$pathway_summary
    background_summary <- background_pathways$pathway_summary
    
    # Merge candidate and background data
    all_pathways_comp <- merge(
      candidate_summary[, c("kegg_id", "pathway_name", "frequency", "unique_loci")],
      background_summary[, c("kegg_id", "frequency", "unique_loci")],
      by = "kegg_id", all = TRUE, suffixes = c("_candidate", "_background")
    )
    
    # Replace NAs with 0
    all_pathways_comp[is.na(all_pathways_comp)] <- 0
    
    # Calculate enrichment metrics
    all_pathways_comp$fold_change <- ifelse(
      all_pathways_comp$frequency_background > 0,
      all_pathways_comp$frequency_candidate / all_pathways_comp$frequency_background,
      ifelse(all_pathways_comp$frequency_candidate > 0, Inf, 1)
    )
    
    all_pathways_comp$enrichment_status <- ifelse(
      all_pathways_comp$fold_change > 2, "Enriched",
      ifelse(all_pathways_comp$fold_change < 0.5, "Depleted", "Similar")
    )
    
    # Sort by fold change
    candidate_pathway_profile <- all_pathways_comp[order(all_pathways_comp$fold_change, decreasing = TRUE), ]
    rownames(candidate_pathway_profile) <- NULL
  }
  
  # 8. Functional clustering
  functional_clusters <- list()
  if (!is.null(network_graph) && igraph::vcount(network_graph) > 2) {
    if (verbose) message("  - Identifying functional pathway clusters...")
    
    # Use community detection to find functional clusters
    tryCatch({
      communities <- igraph::cluster_louvain(network_graph)
      
      # Extract cluster information
      cluster_membership <- igraph::membership(communities)
      
      for (cluster_id in unique(cluster_membership)) {
        cluster_pathways <- names(cluster_membership)[cluster_membership == cluster_id]
        
        if (length(cluster_pathways) >= 2) {  # Only include clusters with multiple pathways
          cluster_info <- pathway_hubs[pathway_hubs$pathway_id %in% cluster_pathways, ]
          
          functional_clusters[[paste0("Cluster_", cluster_id)]] <- list(
            pathway_ids = cluster_pathways,
            pathway_count = length(cluster_pathways),
            pathways = cluster_info,
            total_loci = sum(cluster_info$unique_loci, na.rm = TRUE),
            avg_connectivity = mean(cluster_info$degree, na.rm = TRUE)
          )
        }
      }
    }, error = function(e) {
      if (verbose) message("    - Community detection failed: ", e$message)
    })
  }
  
  # Compile results
  result <- list(
    network_graph = network_graph,
    pathway_hubs = pathway_hubs,
    crosstalk_analysis = crosstalk_analysis,
    candidate_pathway_profile = candidate_pathway_profile,
    network_metrics = network_metrics,
    functional_clusters = functional_clusters,
    analysis_parameters = list(
      candidate_loci_count = nrow(candidate_loci),
      background_loci_count = if (!is.null(background_loci)) nrow(background_loci) else 0,
      include_pathway_crosstalk = include_pathway_crosstalk,
      include_centrality_analysis = include_centrality_analysis,
      min_pathway_size = min_pathway_size,
      max_pathways = max_pathways
    )
  )
  
  if (verbose) {
    message("Pathway network analysis complete:")
    if (!is.null(network_graph)) {
      message("  - Network nodes (pathways): ", igraph::vcount(network_graph))
      message("  - Network edges (interactions): ", igraph::ecount(network_graph))
      message("  - Pathway hubs identified: ", nrow(pathway_hubs))
      message("  - Functional clusters: ", length(functional_clusters))
    } else {
      message("  - No pathway network could be constructed")
    }
    
    if (!is.null(background_loci) && nrow(candidate_pathway_profile) > 0) {
      enriched_count <- sum(candidate_pathway_profile$enrichment_status == "Enriched", na.rm = TRUE)
      message("  - Enriched pathways vs background: ", enriched_count)
    }
  }
  
  return(result)
}