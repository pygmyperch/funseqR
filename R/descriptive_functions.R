#' Descriptive Functions for Functional Annotation Analysis
#'
#' This file contains functions to summarize and describe functional annotations
#' from BLAST results, GO terms, UniProt data, and enrichment analyses.

#' Summarize GO annotations in the database
#'
#' Provides comprehensive statistics about GO term annotations across different
#' ontology categories (BP, MF, CC) and blast parameter sets.
#'
#' @param con Database connection object
#' @param blast_param_id Integer. Specific BLAST parameter set to analyze. If NULL, analyzes all.
#' @param include_evidence Logical. Include GO evidence code statistics. Default is TRUE.
#' @param min_frequency Integer. Minimum frequency threshold for GO terms to include. Default is 1.
#' @param verbose Logical. Print progress information. Default is TRUE.
#'
#' @return List containing GO annotation summaries:
#' \itemize{
#'   \item overview: Summary statistics by ontology category
#'   \item top_terms: Most frequent GO terms by category
#'   \item evidence_summary: GO evidence code distribution (if include_evidence = TRUE)
#'   \item coverage_stats: Coverage statistics across BLAST parameters
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
#' # Analyze specific BLAST parameter set
#' go_summary_blast1 <- summarize_go_annotations(con, blast_param_id = 1)
#' 
#' close_funseq_db(con)
#' }
#'
#' @export
summarize_go_annotations <- function(con, blast_param_id = NULL, include_evidence = TRUE, 
                                   min_frequency = 1, verbose = TRUE) {
  
  if (verbose) message("Summarizing GO annotations...")
  
  # Check if required tables exist
  tables <- DBI::dbListTables(con)
  required_tables <- c("go_terms", "annotations", "blast_results", "blast_parameters")
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
    ", base_where, "
    GROUP BY gt.go_category
    ORDER BY gt.go_category
  ")
  
  overview <- DBI::dbGetQuery(con, overview_query, params)
  
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
      ", base_where, "
      GROUP BY gt.go_category, gt.go_evidence
      ORDER BY gt.go_category, count DESC
    ")
    
    evidence_summary <- DBI::dbGetQuery(con, evidence_query, params)
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
  
  # Compile results
  result <- list(
    overview = overview,
    top_terms = top_terms,
    evidence_summary = evidence_summary,
    coverage_stats = coverage_stats,
    parameters = list(
      blast_param_id = blast_param_id,
      include_evidence = include_evidence,
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
  if (length(base_conditions) > 0 || loci_condition != "") {
    conditions <- c(base_conditions)
    if (loci_condition != "") {
      where_clause <- paste("WHERE", paste(conditions, collapse = " AND "), loci_condition)
    } else {
      where_clause <- paste("WHERE", paste(conditions, collapse = " AND "))
    }
  } else if (loci_condition != "") {
    where_clause <- paste("WHERE", substr(loci_condition, 5, nchar(loci_condition)))  # Remove "AND "
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
  
  loci_summary <- DBI::dbGetQuery(con, loci_summary_query, params)
  
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
    
    pathway_profile <- DBI::dbGetQuery(con, pathway_query, params)
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
    
    blast_quality <- DBI::dbGetQuery(con, blast_quality_query, params)
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
  
  functional_diversity <- DBI::dbGetQuery(con, diversity_query, params)
  
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