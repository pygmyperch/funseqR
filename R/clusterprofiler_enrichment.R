#' clusterProfiler Integration for GO Enrichment Analysis
#'
#' Functions to integrate clusterProfiler's robust enrichment analysis with funseqR workflows

#' Run GO enrichment analysis using clusterProfiler
#'
#' @param con Database connection object
#' @param candidate_file_id Integer. File ID of candidate dataset (e.g., adaptive loci)
#' @param background_file_id Integer. File ID of background dataset (e.g., all sequenced loci)
#' @param blast_param_id Integer. BLAST parameter ID to use for annotations
#' @param ontologies Character vector. GO ontologies to test: c("BP", "MF", "CC"). Default is c("BP", "MF", "CC")
#' @param pvalue_cutoff Numeric. P-value cutoff for significance. Default is 0.05
#' @param padjust_method Character. P-value adjustment method. Default is "BH" (Benjamini-Hochberg)
#' @param min_gs_size Integer. Minimum gene set size for testing. Default is 5
#' @param max_gs_size Integer. Maximum gene set size for testing. Default is 500
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return List containing clusterProfiler enrichResult objects for each ontology
#'
#' @details
#' This function integrates clusterProfiler's robust enrichment analysis with the funseqR workflow.
#' It uses the hypergeometric test (same as the original funseqR implementation) but leverages
#' clusterProfiler's well-tested statistical framework and enhanced visualization capabilities.
#' 
#' The function:
#' 1. Links candidate loci to background annotations by genomic position
#' 2. Extracts UniProt accessions for candidate and background gene sets
#' 3. Prepares GO term annotations in clusterProfiler format
#' 4. Runs enrichment analysis using clusterProfiler's enricher() function
#' 5. Returns results with built-in visualization and summary capabilities
#'
#' @examples
#' \dontrun{
#' # After completing funseqR workflow (VCF import, BLAST, annotation):
#' enrichment_results <- run_clusterprofiler_enrichment(
#'   con = con,
#'   candidate_file_id = vcf_cand_import$file_id,
#'   background_file_id = vcf_import$file_id,
#'   blast_param_id = blast_results$blast_param_id,
#'   ontologies = c("BP", "MF", "CC"),
#'   pvalue_cutoff = 0.1  # More lenient threshold
#' )
#' 
#' # View results
#' summary(enrichment_results$BP)
#' head(enrichment_results$BP@result)
#' 
#' # Create visualizations
#' dotplot(enrichment_results$BP)
#' barplot(enrichment_results$BP)
#' }
#'
#' @export
run_clusterprofiler_enrichment <- function(con, candidate_file_id, background_file_id, blast_param_id,
                                         ontologies = c("BP", "MF", "CC"),
                                         pvalue_cutoff = 0.05,
                                         padjust_method = "BH",
                                         min_gs_size = 5,
                                         max_gs_size = 500,
                                         verbose = TRUE) {
  
  if (verbose) message("=== Starting clusterProfiler GO Enrichment Analysis ===")
  
  # Check required packages
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("clusterProfiler package is required. Install with: BiocManager::install('clusterProfiler')")
  }
  if (!requireNamespace("DOSE", quietly = TRUE)) {
    stop("DOSE package is required. Install with: BiocManager::install('DOSE')")
  }
  
  # Step 1: Extract candidate gene list
  if (verbose) message("Step 1: Extracting candidate gene list...")
  candidate_genes <- .extract_candidate_uniprot_ids(con, candidate_file_id, background_file_id, blast_param_id, verbose)
  
  if (length(candidate_genes) == 0) {
    stop("No candidate genes found. Check that candidate and background files are properly linked.")
  }
  
  # Step 2: Extract background gene list
  if (verbose) message("Step 2: Extracting background gene list...")
  background_genes <- .extract_background_uniprot_ids(con, background_file_id, blast_param_id, verbose)
  
  if (length(background_genes) == 0) {
    stop("No background genes found. Check that background file has GO annotations.")
  }
  
  if (verbose) {
    message("  - Candidate genes: ", length(candidate_genes))
    message("  - Background genes: ", length(background_genes))
    message("  - Overlap: ", length(intersect(candidate_genes, background_genes)))
  }
  
  # Step 3: Prepare GO annotation data
  if (verbose) message("Step 3: Preparing GO annotation data...")
  go_data <- .prepare_clusterprofiler_go_data(con, background_file_id, blast_param_id, verbose)
  
  # Step 4: Run enrichment analysis for each ontology
  if (verbose) message("Step 4: Running enrichment analysis...")
  results <- list()
  
  for (ontology in ontologies) {
    if (verbose) message("  - Analyzing ", ontology, " ontology...")
    
    # Filter GO data for this ontology
    ontology_data <- .filter_go_data_by_ontology(go_data, ontology)
    
    if (nrow(ontology_data$term2gene) == 0) {
      if (verbose) message("    - No ", ontology, " terms found, skipping...")
      results[[ontology]] <- NULL
      next
    }
    
    # Run clusterProfiler enrichment
    tryCatch({
      enrichment_result <- clusterProfiler::enricher(
        gene = candidate_genes,
        universe = background_genes,
        TERM2GENE = ontology_data$term2gene,
        TERM2NAME = ontology_data$term2name,
        pvalueCutoff = pvalue_cutoff,
        pAdjustMethod = padjust_method,
        minGSSize = min_gs_size,
        maxGSSize = max_gs_size
      )
      
      results[[ontology]] <- enrichment_result
      
      if (verbose) {
        if (is.null(enrichment_result) || nrow(enrichment_result@result) == 0) {
          message("    - No significant terms found")
        } else {
          significant_count <- sum(enrichment_result@result$p.adjust < pvalue_cutoff)
          message("    - Terms tested: ", nrow(enrichment_result@result))
          message("    - Significant terms: ", significant_count)
        }
      }
      
    }, error = function(e) {
      warning("Error in ", ontology, " enrichment: ", e$message)
      results[[ontology]] <- NULL
    })
  }
  
  # Step 5: Create summary
  if (verbose) message("Step 5: Creating analysis summary...")
  summary_info <- .create_enrichment_summary(results, candidate_genes, background_genes, 
                                            candidate_file_id, background_file_id, blast_param_id)
  
  if (verbose) {
    message("=== Analysis Complete ===")
    message("Candidate genes analyzed: ", length(candidate_genes))
    message("Background genes: ", length(background_genes))
    total_significant <- sum(sapply(results, function(x) {
      if (is.null(x)) return(0)
      sum(x@result$p.adjust < pvalue_cutoff, na.rm = TRUE)
    }))
    message("Total significant terms: ", total_significant)
  }
  
  # Add summary to results
  results$summary <- summary_info
  results$parameters <- list(
    candidate_file_id = candidate_file_id,
    background_file_id = background_file_id,
    blast_param_id = blast_param_id,
    ontologies = ontologies,
    pvalue_cutoff = pvalue_cutoff,
    padjust_method = padjust_method,
    min_gs_size = min_gs_size,
    max_gs_size = max_gs_size
  )
  
  class(results) <- c("funseqR_enrichment", "list")
  return(results)
}

#' Extract candidate UniProt IDs by linking positions
#' @keywords internal
.extract_candidate_uniprot_ids <- function(con, candidate_file_id, background_file_id, blast_param_id, verbose = FALSE) {
  
  if (verbose) message("  - Linking candidate positions to background annotations...")
  
  # Query to link candidate positions to background annotations
  # This replicates the logic from the old workflow
  query <- "
    SELECT DISTINCT a.uniprot_accession
    FROM vcf_data c
    JOIN vcf_data r ON (c.chromosome = r.chromosome AND c.position = r.position)
    JOIN flanking_sequences fs ON r.vcf_id = fs.vcf_id
    JOIN blast_results br ON fs.flanking_id = br.flanking_id
    JOIN annotations a ON br.blast_result_id = a.blast_result_id
    WHERE c.file_id = ? AND r.file_id = ? AND br.blast_param_id = ?
      AND a.uniprot_accession IS NOT NULL
      AND a.uniprot_accession != ''
  "
  
  result <- DBI::dbGetQuery(con, query, list(candidate_file_id, background_file_id, blast_param_id))
  
  if (verbose) message("    - Found ", nrow(result), " candidate genes")
  
  return(result$uniprot_accession)
}

#' Extract background UniProt IDs with GO annotations
#' @keywords internal
.extract_background_uniprot_ids <- function(con, background_file_id, blast_param_id, verbose = FALSE) {
  
  if (verbose) message("  - Extracting background genes with GO annotations...")
  
  # Query to get all background genes that have GO annotations
  query <- "
    SELECT DISTINCT a.uniprot_accession
    FROM vcf_data v
    JOIN flanking_sequences fs ON v.vcf_id = fs.vcf_id
    JOIN blast_results br ON fs.flanking_id = br.flanking_id
    JOIN annotations a ON br.blast_result_id = a.blast_result_id
    JOIN go_terms gt ON a.annotation_id = gt.annotation_id
    WHERE v.file_id = ? AND br.blast_param_id = ?
      AND a.uniprot_accession IS NOT NULL
      AND a.uniprot_accession != ''
  "
  
  result <- DBI::dbGetQuery(con, query, list(background_file_id, blast_param_id))
  
  if (verbose) message("    - Found ", nrow(result), " background genes")
  
  return(result$uniprot_accession)
}

#' Prepare GO annotation data for clusterProfiler
#' @keywords internal
.prepare_clusterprofiler_go_data <- function(con, background_file_id, blast_param_id, verbose = FALSE) {
  
  if (verbose) message("  - Preparing GO term mappings...")
  
  # Query to get gene-to-GO and GO-to-name mappings
  query <- "
    SELECT DISTINCT 
      gt.go_id,
      gt.go_term,
      gt.go_category,
      a.uniprot_accession
    FROM vcf_data v
    JOIN flanking_sequences fs ON v.vcf_id = fs.vcf_id
    JOIN blast_results br ON fs.flanking_id = br.flanking_id
    JOIN annotations a ON br.blast_result_id = a.blast_result_id
    JOIN go_terms gt ON a.annotation_id = gt.annotation_id
    WHERE v.file_id = ? AND br.blast_param_id = ?
      AND a.uniprot_accession IS NOT NULL
      AND a.uniprot_accession != ''
      AND gt.go_id IS NOT NULL
      AND gt.go_id != ''
  "
  
  go_annotations <- DBI::dbGetQuery(con, query, list(background_file_id, blast_param_id))
  
  if (nrow(go_annotations) == 0) {
    stop("No GO annotations found for background dataset")
  }
  
  # Create TERM2GENE mapping (GO ID -> Gene)
  term2gene <- go_annotations[, c("go_id", "uniprot_accession")]
  colnames(term2gene) <- c("term", "gene")
  
  # Create TERM2NAME mapping (GO ID -> GO term name)
  term2name <- unique(go_annotations[, c("go_id", "go_term")])
  colnames(term2name) <- c("term", "name")
  
  if (verbose) {
    message("    - Total GO terms: ", length(unique(go_annotations$go_id)))
    message("    - Total gene-term associations: ", nrow(term2gene))
  }
  
  return(list(
    term2gene = term2gene,
    term2name = term2name,
    raw_data = go_annotations
  ))
}

#' Filter GO data by ontology
#' @keywords internal
.filter_go_data_by_ontology <- function(go_data, ontology) {
  
  # Map ontology codes
  ontology_map <- c("BP" = "P", "MF" = "F", "CC" = "C")
  category_code <- ontology_map[ontology]
  
  if (is.na(category_code)) {
    stop("Invalid ontology. Must be 'BP', 'MF', or 'CC'")
  }
  
  # Filter raw data by category
  ontology_annotations <- go_data$raw_data[go_data$raw_data$go_category == category_code, ]
  
  if (nrow(ontology_annotations) == 0) {
    return(list(term2gene = data.frame(term = character(0), gene = character(0)),
                term2name = data.frame(term = character(0), name = character(0))))
  }
  
  # Create filtered mappings
  term2gene <- ontology_annotations[, c("go_id", "uniprot_accession")]
  colnames(term2gene) <- c("term", "gene")
  
  term2name <- unique(ontology_annotations[, c("go_id", "go_term")])
  colnames(term2name) <- c("term", "name")
  
  return(list(
    term2gene = term2gene,
    term2name = term2name
  ))
}

#' Create enrichment analysis summary
#' @keywords internal
.create_enrichment_summary <- function(results, candidate_genes, background_genes, 
                                     candidate_file_id, background_file_id, blast_param_id) {
  
  summary_stats <- data.frame(
    Ontology = names(results)[names(results) != "summary"],
    stringsAsFactors = FALSE
  )
  
  summary_stats$Terms_Tested <- sapply(summary_stats$Ontology, function(ont) {
    if (is.null(results[[ont]])) return(0)
    nrow(results[[ont]]@result)
  })
  
  summary_stats$Significant_Terms <- sapply(summary_stats$Ontology, function(ont) {
    if (is.null(results[[ont]])) return(0)
    sum(results[[ont]]@result$p.adjust < 0.05, na.rm = TRUE)
  })
  
  summary_stats$Top_Enrichment <- sapply(summary_stats$Ontology, function(ont) {
    if (is.null(results[[ont]]) || nrow(results[[ont]]@result) == 0) return(0)
    # Calculate fold enrichment (Count/Expected ratio)
    result_df <- results[[ont]]@result
    if ("Count" %in% colnames(result_df) && "Expected" %in% colnames(result_df)) {
      max_enrich <- max(result_df$Count / result_df$Expected, na.rm = TRUE)
      return(round(max_enrich, 2))
    }
    return(0)
  })
  
  return(list(
    analysis_date = Sys.time(),
    candidate_file_id = candidate_file_id,
    background_file_id = background_file_id,
    blast_param_id = blast_param_id,
    candidate_genes = length(candidate_genes),
    background_genes = length(background_genes),
    summary_stats = summary_stats
  ))
}

#' Print method for funseqR enrichment results
#' @export
print.funseqR_enrichment <- function(x, ...) {
  cat("=== funseqR clusterProfiler Enrichment Results ===\n\n")
  
  if (!is.null(x$summary)) {
    cat("Analysis Date:", format(x$summary$analysis_date), "\n")
    cat("Candidate Genes:", x$summary$candidate_genes, "\n")
    cat("Background Genes:", x$summary$background_genes, "\n\n")
    
    if (!is.null(x$summary$summary_stats)) {
      cat("Results Summary:\n")
      print(x$summary$summary_stats)
      cat("\n")
    }
  }
  
  # Show available result objects
  result_objects <- names(x)[!names(x) %in% c("summary", "parameters")]
  if (length(result_objects) > 0) {
    cat("Available result objects:", paste(result_objects, collapse = ", "), "\n")
    cat("Access results with: $BP, $MF, $CC\n")
    cat("Create plots with: dotplot(results$BP), barplot(results$BP)\n")
  }
  
  invisible(x)
}