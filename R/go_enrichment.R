#' GO Enrichment Analysis Functions
#'
#' This file contains functions for performing Gene Ontology (GO) enrichment analysis
#' on candidate loci vs background datasets, with visualization capabilities.
#'

#' Import candidate adaptive loci and link to existing annotations
#'
#' @param con Database connection object
#' @param project_id Integer. Project ID in the database
#' @param candidate_vcf_file Character. Path to candidate VCF file
#' @param background_file_id Integer. File ID of the background/reference dataset
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return List containing file_id, bed_file path, and annotation linkage results
#'
#' @details
#' This function imports a candidate VCF file and links the variants to existing
#' annotations in the database by matching genomic positions. It's designed for
#' comparative analysis between candidate adaptive loci and a larger background dataset.
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#' candidate_import <- import_candidate_loci(con, 1, "candidates.vcf", 1)
#' }
#'
#' @export
import_candidate_loci <- function(con, project_id, candidate_vcf_file, background_file_id, verbose = TRUE) {
  
  if (verbose) message("Importing candidate loci from: ", candidate_vcf_file)
  
  # Import candidate VCF using existing function
  candidate_import <- import_vcf_to_db(con, project_id, candidate_vcf_file)
  
  if (verbose) message("Creating BED file for candidate loci...")
  
  # Create BED file for candidate loci
  candidate_bed <- vcf2bed_db(con, candidate_import$file_id)
  
  if (verbose) message("Linking candidates to existing annotations...")
  
  # Link candidates to existing annotations via genomic overlap
  candidate_annotations <- link_candidates_to_annotations(con, candidate_import$file_id, background_file_id, verbose = verbose)
  
  if (verbose) {
    message("Import complete:")
    message("  - Candidate file ID: ", candidate_import$file_id)
    message("  - BED file created: ", candidate_bed)
    message("  - Linked annotations: ", nrow(candidate_annotations))
  }
  
  return(list(
    file_id = candidate_import$file_id,
    bed_file = candidate_bed,
    linked_annotations = candidate_annotations,
    import_summary = candidate_import
  ))
}

#' Link candidate loci to existing annotations by genomic position
#'
#' @param con Database connection object
#' @param candidate_file_id Integer. File ID of candidate dataset
#' @param background_file_id Integer. File ID of background dataset
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return Data frame with candidate variants linked to annotations
#'
#' @details
#' Links candidate variants to existing functional annotations by matching
#' chromosome and position coordinates. This assumes that both datasets
#' use the same reference genome and coordinate system.
#'
#' @export
link_candidates_to_annotations <- function(con, candidate_file_id, background_file_id, verbose = TRUE) {
  
  if (verbose) message("Querying database for annotation linkages...")
  
  query <- "
    SELECT DISTINCT 
      c.vcf_id as candidate_vcf_id,
      c.chromosome,
      c.position,
      c.ref as candidate_ref,
      c.alt as candidate_alt,
      a.annotation_id,
      a.uniprot_accession,
      a.gene_names,
      br.e_value,
      br.bit_score,
      br.percent_identity
    FROM vcf_data c
    JOIN vcf_data r ON (c.chromosome = r.chromosome 
                       AND c.position = r.position)
    JOIN flanking_sequences fs ON r.vcf_id = fs.vcf_id
    JOIN blast_results br ON fs.flanking_id = br.flanking_id  
    JOIN annotations a ON br.blast_result_id = a.blast_result_id
    WHERE c.file_id = ? AND r.file_id = ?
    ORDER BY c.chromosome, c.position
  "
  
  result <- DBI::dbGetQuery(con, query, params = list(candidate_file_id, background_file_id))
  
  if (verbose) {
    message("Linkage complete:")
    message("  - Candidate variants with annotations: ", length(unique(result$candidate_vcf_id)))
    message("  - Total annotation links: ", nrow(result))
    message("  - Unique proteins: ", length(unique(result$uniprot_accession)))
  }
  
  return(result)
}

#' Extract GO terms for foreground and background gene sets
#'
#' @param con Database connection object
#' @param foreground_file_id Integer. File ID of candidate/foreground dataset
#' @param background_file_id Integer. File ID of background dataset
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return List containing foreground and background GO term data
#'
#' @details
#' Extracts GO terms associated with proteins from both foreground (candidate)
#' and background datasets. Creates gene-to-GO mappings required for enrichment testing.
#'
#' @examples
#' \dontrun{
#' go_data <- extract_go_terms_for_enrichment(con, candidate_file_id, background_file_id)
#' str(go_data)
#' }
#'
#' @export
extract_go_terms_for_enrichment <- function(con, foreground_file_id, background_file_id, verbose = TRUE) {
  
  if (verbose) message("Extracting GO terms for enrichment analysis...")
  
  # First check if foreground file has direct annotations (complete analysis)
  # or if it's a candidate file that needs to use linkage data
  foreground_direct_count <- DBI::dbGetQuery(con, "
    SELECT COUNT(*) as count FROM vcf_data v
    JOIN flanking_sequences fs ON v.vcf_id = fs.vcf_id
    JOIN blast_results br ON fs.flanking_id = br.flanking_id
    JOIN annotations a ON br.blast_result_id = a.blast_result_id
    WHERE v.file_id = ?
  ", params = list(foreground_file_id))$count
  
  if (foreground_direct_count > 0) {
    # Standard query for datasets with their own annotations
    go_query <- "
      SELECT DISTINCT 
        a.uniprot_accession,
        gt.go_id,
        gt.go_term,
        gt.go_category,
        gt.go_evidence
      FROM vcf_data v
      JOIN flanking_sequences fs ON v.vcf_id = fs.vcf_id
      JOIN blast_results br ON fs.flanking_id = br.flanking_id
      JOIN annotations a ON br.blast_result_id = a.blast_result_id
      JOIN go_terms gt ON a.annotation_id = gt.annotation_id
      WHERE v.file_id = ?
    "
    
    if (verbose) message("  - Extracting foreground GO terms...")
    foreground_go <- DBI::dbGetQuery(con, go_query, params = list(foreground_file_id))
  } else {
    # For candidate files, use the linked annotations
    if (verbose) message("  - Extracting foreground GO terms via linkage...")
    
    foreground_go_query <- "
      SELECT DISTINCT 
        a.uniprot_accession,
        gt.go_id,
        gt.go_term,
        gt.go_category,
        gt.go_evidence
      FROM vcf_data c
      JOIN vcf_data r ON (c.chromosome = r.chromosome 
                         AND c.position = r.position)
      JOIN flanking_sequences fs ON r.vcf_id = fs.vcf_id
      JOIN blast_results br ON fs.flanking_id = br.flanking_id  
      JOIN annotations a ON br.blast_result_id = a.blast_result_id
      JOIN go_terms gt ON a.annotation_id = gt.annotation_id
      WHERE c.file_id = ? AND r.file_id = ?
    "
    
    foreground_go <- DBI::dbGetQuery(con, foreground_go_query, 
                                   params = list(foreground_file_id, background_file_id))
  }
  
  # Background always uses standard query
  background_go_query <- "
    SELECT DISTINCT 
      a.uniprot_accession,
      gt.go_id,
      gt.go_term,
      gt.go_category,
      gt.go_evidence
    FROM vcf_data v
    JOIN flanking_sequences fs ON v.vcf_id = fs.vcf_id
    JOIN blast_results br ON fs.flanking_id = br.flanking_id
    JOIN annotations a ON br.blast_result_id = a.blast_result_id
    JOIN go_terms gt ON a.annotation_id = gt.annotation_id
    WHERE v.file_id = ?
  "
  
  if (verbose) message("  - Extracting background GO terms...")
  background_go <- DBI::dbGetQuery(con, background_go_query, params = list(background_file_id))
  
  # Create gene-to-GO mapping lists
  foreground_gene2go <- split(foreground_go$go_id, foreground_go$uniprot_accession)
  background_gene2go <- split(background_go$go_id, background_go$uniprot_accession)
  
  # Get unique GO terms and their details
  all_go_terms <- rbind(
    foreground_go[, c("go_id", "go_term", "go_category")],
    background_go[, c("go_id", "go_term", "go_category")]
  )
  all_go_terms <- all_go_terms[!duplicated(all_go_terms), ]
  
  if (verbose) {
    message("GO term extraction complete:")
    message("  - Foreground genes: ", length(unique(foreground_go$uniprot_accession)))
    message("  - Background genes: ", length(unique(background_go$uniprot_accession)))
    message("  - Total unique GO terms: ", nrow(all_go_terms))
    message("  - Biological Process terms: ", sum(all_go_terms$go_category == "P"))
    message("  - Molecular Function terms: ", sum(all_go_terms$go_category == "F"))
    message("  - Cellular Component terms: ", sum(all_go_terms$go_category == "C"))
  }
  
  return(list(
    foreground = list(
      genes = unique(foreground_go$uniprot_accession),
      gene2go = foreground_gene2go,
      go_terms = foreground_go
    ),
    background = list(
      genes = unique(background_go$uniprot_accession),
      gene2go = background_gene2go,
      go_terms = background_go
    ),
    all_go_terms = all_go_terms
  ))
}

#' Perform GO enrichment analysis using hypergeometric test
#'
#' @param go_data List. Output from extract_go_terms_for_enrichment()
#' @param ontology Character. GO ontology: "BP" (Biological Process), "MF" (Molecular Function), or "CC" (Cellular Component)
#' @param min_genes Integer. Minimum genes required for a GO term to be tested. Default is 5
#' @param max_genes Integer. Maximum genes for a GO term (to exclude very broad terms). Default is 500
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return Data frame with enrichment results, sorted by adjusted p-value
#'
#' @details
#' Performs hypergeometric enrichment testing for GO terms. Tests whether each GO term
#' is overrepresented in the foreground set compared to the background set.
#' Applies FDR correction for multiple testing.
#'
#' @examples
#' \dontrun{
#' go_data <- extract_go_terms_for_enrichment(con, fg_id, bg_id)
#' bp_results <- perform_go_enrichment(go_data, "BP")
#' head(bp_results)
#' }
#'
#' @export
perform_go_enrichment <- function(go_data, ontology = "BP", min_genes = 5, max_genes = 500, verbose = TRUE) {
  
  if (verbose) message("Performing GO enrichment analysis for ontology: ", ontology)
  
  # Map ontology codes
  ontology_map <- c("BP" = "P", "MF" = "F", "CC" = "C")
  category_code <- ontology_map[ontology]
  
  if (is.na(category_code)) {
    stop("Invalid ontology. Must be 'BP', 'MF', or 'CC'")
  }
  
  # Filter GO terms by category
  relevant_terms <- go_data$all_go_terms[go_data$all_go_terms$go_category == category_code, "go_id"]
  relevant_terms <- unique(relevant_terms)
  
  if (verbose) message("  - Testing ", length(relevant_terms), " GO terms in category ", ontology)
  
  # Calculate total gene counts
  total_fg <- length(go_data$foreground$genes)
  total_bg <- length(go_data$background$genes)
  
  if (verbose) message("  - Foreground genes: ", total_fg, ", Background genes: ", total_bg)
  
  # Calculate enrichment for each GO term
  enrichment_list <- list()
  
  for (go_term in relevant_terms) {
    
    # Count genes with this GO term in each set
    fg_with_term <- sum(sapply(go_data$foreground$gene2go, function(x) go_term %in% x))
    bg_with_term <- sum(sapply(go_data$background$gene2go, function(x) go_term %in% x))
    
    # Apply gene count filters
    if (bg_with_term < min_genes || bg_with_term > max_genes) {
      next
    }
    
    # Skip if no foreground genes have this term
    if (fg_with_term == 0) {
      next
    }
    
    # Hypergeometric test
    # P(X >= fg_with_term) where X ~ Hypergeometric(total_bg, bg_with_term, total_fg)
    p_value <- phyper(fg_with_term - 1, bg_with_term, total_bg - bg_with_term, total_fg, lower.tail = FALSE)
    
    # Calculate expected count and fold enrichment
    expected <- (bg_with_term / total_bg) * total_fg
    fold_enrichment <- ifelse(expected > 0, fg_with_term / expected, Inf)
    
    # Get GO term details
    term_info <- go_data$all_go_terms[go_data$all_go_terms$go_id == go_term, ][1, ]
    
    enrichment_list[[length(enrichment_list) + 1]] <- data.frame(
      go_id = go_term,
      go_term = term_info$go_term,
      go_category = ontology,
      foreground_count = fg_with_term,
      background_count = bg_with_term,
      total_foreground = total_fg,
      total_background = total_bg,
      expected_count = expected,
      fold_enrichment = fold_enrichment,
      p_value = p_value,
      stringsAsFactors = FALSE
    )
  }
  
  # Combine all results
  if (length(enrichment_list) > 0) {
    enrichment_results <- do.call(rbind, enrichment_list)
  } else {
    enrichment_results <- data.frame()
  }
  
  # Handle case where no terms pass filters
  if (is.null(enrichment_results) || nrow(enrichment_results) == 0) {
    if (verbose) message("  - No GO terms passed filtering criteria")
    return(data.frame())
  }
  
  # Multiple testing correction
  enrichment_results$p_adjusted <- p.adjust(enrichment_results$p_value, method = "fdr")
  
  # Add significance levels
  enrichment_results$significance_level <- ifelse(
    enrichment_results$p_adjusted < 0.01, "highly_significant",
    ifelse(enrichment_results$p_adjusted < 0.05, "significant",
           ifelse(enrichment_results$p_adjusted < 0.1, "trending", "not_significant"))
  )
  
  # Sort by significance and fold enrichment
  enrichment_results <- enrichment_results[order(enrichment_results$p_adjusted, -enrichment_results$fold_enrichment), ]
  
  if (verbose) {
    sig_count <- sum(enrichment_results$p_adjusted < 0.05)
    message("  - Enrichment analysis complete: ", nrow(enrichment_results), " terms tested, ", sig_count, " significantly enriched")
  }
  
  return(enrichment_results)
}

#' Store GO enrichment analysis results in the database
#'
#' @param con Database connection object
#' @param project_id Integer. Project ID
#' @param foreground_file_id Integer. File ID of foreground dataset
#' @param background_file_id Integer. File ID of background dataset
#' @param enrichment_results Data frame. Results from perform_go_enrichment()
#' @param ontology Character. GO ontology tested
#' @param parameters List. Analysis parameters for reproducibility
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return Integer. The enrichment_id of the stored analysis
#'
#' @export
store_go_enrichment_results <- function(con, project_id, foreground_file_id, background_file_id, 
                                       enrichment_results, ontology, parameters = NULL, verbose = TRUE) {
  
  if (verbose) message("Storing GO enrichment results in database...")
  
  # Ensure tables exist
  tables <- DBI::dbListTables(con)
  if (!all(c("go_enrichment_analyses", "go_enrichment_results") %in% tables)) {
    stop("GO enrichment tables not found in database. Please upgrade schema.")
  }
  
  # Prepare analysis parameters
  if (is.null(parameters)) {
    parameters <- list()
  }
  params_json <- jsonlite::toJSON(parameters, auto_unbox = TRUE)
  
  # Insert analysis record
  analysis_query <- "
    INSERT INTO go_enrichment_analyses 
    (project_id, foreground_file_id, background_file_id, ontology, analysis_date, 
     total_foreground_genes, total_background_genes, analysis_parameters)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
  "
  
  DBI::dbExecute(con, analysis_query, params = list(
    project_id,
    foreground_file_id,
    background_file_id,
    ontology,
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    if(nrow(enrichment_results) > 0) enrichment_results$total_foreground[1] else 0,
    if(nrow(enrichment_results) > 0) enrichment_results$total_background[1] else 0,
    params_json
  ))
  
  # Get the analysis ID
  enrichment_id <- DBI::dbGetQuery(con, "SELECT last_insert_rowid() as id")$id
  
  # Insert results if any
  if (nrow(enrichment_results) > 0) {
    
    # Prepare results for insertion
    results_to_insert <- enrichment_results
    results_to_insert$enrichment_id <- enrichment_id
    
    # Select only the columns we need in the right order
    results_to_insert <- results_to_insert[, c("enrichment_id", "go_id", "go_term", "go_category", 
                                             "foreground_count", "background_count", "total_foreground", 
                                             "total_background", "expected_count", "fold_enrichment", 
                                             "p_value", "p_adjusted", "significance_level")]
    
    # Insert results row by row (more reliable than batch insert)
    result_query <- "
      INSERT INTO go_enrichment_results 
      (enrichment_id, go_id, go_term, go_category, foreground_count, background_count,
       total_foreground, total_background, expected_count, fold_enrichment, 
       p_value, p_adjusted, significance_level)
      VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    "
    
    for (i in 1:nrow(results_to_insert)) {
      row_data <- as.list(results_to_insert[i, ])
      DBI::dbExecute(con, result_query, params = row_data)
    }
  }
  
  if (verbose) {
    message("Storage complete:")
    message("  - Analysis ID: ", enrichment_id)
    message("  - Results stored: ", nrow(enrichment_results))
  }
  
  return(enrichment_id)
}

#' Retrieve stored GO enrichment results
#'
#' @param con Database connection object
#' @param enrichment_id Integer. Analysis ID to retrieve
#' @param significance_filter Character. Filter by significance level. Default is NULL (no filter)
#'
#' @return List containing analysis metadata and results
#'
#' @export
get_go_enrichment_results <- function(con, enrichment_id, significance_filter = NULL) {
  
  # Get analysis metadata
  analysis_query <- "SELECT * FROM go_enrichment_analyses WHERE enrichment_id = ?"
  analysis_info <- DBI::dbGetQuery(con, analysis_query, params = list(enrichment_id))
  
  if (nrow(analysis_info) == 0) {
    stop("No enrichment analysis found with ID: ", enrichment_id)
  }
  
  # Get results
  results_query <- "
    SELECT * FROM go_enrichment_results 
    WHERE enrichment_id = ?
    ORDER BY p_adjusted, fold_enrichment DESC
  "
  
  results <- DBI::dbGetQuery(con, results_query, params = list(enrichment_id))
  
  # Apply significance filter if requested
  if (!is.null(significance_filter)) {
    results <- results[results$significance_level %in% significance_filter, ]
  }
  
  return(list(
    analysis_info = as.list(analysis_info[1, ]),
    results = results
  ))
}