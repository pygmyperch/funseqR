#' Descriptive Functions for Functional Annotation Analysis
#'
#' This file contains functions to summarize and describe functional annotations
#' from BLAST results, GO terms, UniProt data, and enrichment analyses.

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
#'
#' @param con Database connection object
#' @param blast_param_id Integer. Specific BLAST parameter set to analyze. If NULL, analyzes all.
#' @param candidate_loci Data frame with chromosome and position columns, or NULL for all loci.
#' @param include_modules Logical. Group pathways into functional modules. Default is TRUE.
#' @param min_frequency Integer. Minimum frequency threshold for pathways to include. Default is 1.
#' @param verbose Logical. Print progress information. Default is TRUE.
#'
#' @return List containing KEGG pathway summaries:
#' \itemize{
#'   \item overview: Summary statistics of pathway coverage
#'   \item top_pathways: Most frequent pathways 
#'   \item pathway_modules: Pathways grouped by functional modules (if include_modules = TRUE)
#'   \item coverage_stats: Coverage statistics across BLAST parameters
#' }
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#' 
#' # Analyze all KEGG pathways with module grouping
#' kegg_summary <- summarize_kegg_pathways(con, include_modules = TRUE)
#' print(kegg_summary$pathway_modules)
#' 
#' # Analyze specific candidate loci
#' candidates <- data.frame(chromosome = c("LG1", "LG2"), position = c(12345, 67890))
#' kegg_candidates <- summarize_kegg_pathways(con, candidate_loci = candidates)
#' 
#' close_funseq_db(con)
#' }
#'
#' @export
summarize_kegg_pathways <- function(con, blast_param_id = NULL, candidate_loci = NULL,
                                   include_modules = TRUE, min_frequency = 1, verbose = TRUE) {
  
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
  top_pathways <- DBI::dbGetQuery(con, top_pathways_query, top_pathways_params)
  
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
    parameters = list(
      blast_param_id = blast_param_id,
      candidate_loci = !is.null(candidate_loci),
      include_modules = include_modules,
      min_frequency = min_frequency
    )
  )
  
  if (verbose) {
    message("KEGG pathway summary complete:")
    if (nrow(overview) > 0 && overview$total_pathway_assignments[1] > 0) {
      message("  - Unique pathways: ", overview$unique_pathways[1])
      message("  - Total assignments: ", overview$total_pathway_assignments[1])
      message("  - Annotated loci: ", overview$annotated_loci[1])
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
  
  if (include_go && !is.null(result$go_modules)) {
    overview_stats$go_module_count <- length(result$go_modules)
    overview_stats$go_total_terms <- sum(sapply(result$go_modules, function(m) m$term_count))
  }
  
  if (include_kegg && !is.null(result$kegg_modules)) {
    overview_stats$kegg_module_count <- length(result$kegg_modules)
    overview_stats$kegg_total_pathways <- sum(sapply(result$kegg_modules, function(m) m$pathway_count))
  }
  
  if (include_cog && !is.null(result$cog_modules)) {
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
    verbose = FALSE
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
  
  overview <- DBI::dbGetQuery(con, overview_query, params)
  
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
  top_proteins <- DBI::dbGetQuery(con, top_proteins_query, top_proteins_params)
  
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
    gene_name_stats <- DBI::dbGetQuery(con, gene_stats_query, gene_stats_params)
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
  
  annotation_quality <- DBI::dbGetQuery(con, quality_query, params)
  
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
  
  if (length(params) > 0) {
    loci_summary <- DBI::dbGetQuery(con, loci_summary_query, params)
  } else {
    loci_summary <- DBI::dbGetQuery(con, loci_summary_query)
  }
  
  # 2. GO profile (using existing function)
  if (verbose) message("  - Analyzing GO term profile...")
  go_profile <- summarize_go_annotations(con, blast_param_id = blast_param_id, verbose = FALSE)
  
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
    
    if (length(params) > 0) {
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
    
    if (length(params) > 0) {
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
  
  if (length(params) > 0) {
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
  category_distribution <- DBI::dbGetQuery(con, distribution_query, distribution_params)
  
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
  
  annotation_coverage <- DBI::dbGetQuery(con, coverage_query, coverage_params)
  
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