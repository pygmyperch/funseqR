# EXPORTED

#' GO Enrichment Analysis Functions
#'
#' This file contains functions for performing Gene Ontology (GO) enrichment analysis
#' on candidate loci vs background datasets, with visualization capabilities.
#'


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
#' @keywords internal
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

  result <- DBI::dbGetQuery(con, query, list(candidate_file_id, background_file_id))

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
#' @param blast_param_id Integer. Optional. Specific BLAST run ID to use for both datasets.
#'   If NULL, uses all available annotations. Default is NULL.
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return List containing foreground and background GO term data
#'
#' @details
#' Extracts GO terms associated with proteins from both foreground (candidate)
#' and background datasets. Creates gene-to-GO mappings required for enrichment testing.
#' 
#' When blast_param_id is specified, only annotations from that specific BLAST run
#' are used for both datasets. This ensures methodological consistency and allows
#' comparison of different search strategies (e.g., ORF sequences vs raw sequences).
#'
#' @examples
#' \dontrun{
#' # Use all available annotations
#' go_data <- extract_go_terms_for_enrichment(con, foreground_file_id, background_file_id)
#' 
#' # Use specific BLAST run (e.g., ORF-based search)
#' go_data_orf <- extract_go_terms_for_enrichment(con, foreground_file_id, background_file_id,
#'                                                blast_param_id = 1)
#' str(go_data_orf)
#' }
#'
#' @keywords internal
extract_go_terms_for_enrichment <- function(con, foreground_file_id, background_file_id, 
                                           blast_param_id = NULL, verbose = TRUE) {

  if (verbose) message("Extracting GO terms for enrichment analysis...")

  # Handle stored candidates case
  if (foreground_file_id == "stored") {
    if (verbose) message("  - Using stored candidate loci from database")
    
    # Check if candidates exist
    candidate_count <- DBI::dbGetQuery(con, "SELECT COUNT(*) as count FROM candidate_loci")$count
    if (candidate_count == 0) {
      stop("No candidate loci found in database. Run define_locus_statistics() or define_candidate_loci() first.")
    }
    
    # Extract GO annotations for stored candidates
    result <- .extract_stored_candidate_go_annotations(con, background_file_id, blast_param_id, verbose)
    return(result)
  }

  # First check if foreground file has direct annotations (complete analysis)
  # or if it's a candidate file that needs to use linkage data
  if (is.null(blast_param_id)) {
    foreground_direct_count <- DBI::dbGetQuery(con, "
      SELECT COUNT(*) as count FROM vcf_data v
      JOIN flanking_sequences fs ON v.vcf_id = fs.vcf_id
      JOIN blast_results br ON fs.flanking_id = br.flanking_id
      JOIN annotations a ON br.blast_result_id = a.blast_result_id
      WHERE v.file_id = ?
    ", list(foreground_file_id))$count
  } else {
    foreground_direct_count <- DBI::dbGetQuery(con, "
      SELECT COUNT(*) as count FROM vcf_data v
      JOIN flanking_sequences fs ON v.vcf_id = fs.vcf_id
      JOIN blast_results br ON fs.flanking_id = br.flanking_id
      JOIN annotations a ON br.blast_result_id = a.blast_result_id
      WHERE v.file_id = ? AND br.blast_param_id = ?
    ", list(foreground_file_id, blast_param_id))$count
  }

  if (foreground_direct_count > 0) {
    # Standard query for datasets with their own annotations
    if (is.null(blast_param_id)) {
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
      if (verbose) message("  - Extracting foreground GO terms (all BLAST runs)...")
      foreground_go <- DBI::dbGetQuery(con, go_query, list(foreground_file_id))
    } else {
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
        WHERE v.file_id = ? AND br.blast_param_id = ?
      "
      if (verbose) message("  - Extracting foreground GO terms (BLAST run ID: ", blast_param_id, ")...")
      foreground_go <- DBI::dbGetQuery(con, go_query, list(foreground_file_id, blast_param_id))
    }
  } else {
    # For candidate files, use the linked annotations
    if (is.null(blast_param_id)) {
      if (verbose) message("  - Extracting foreground GO terms via linkage (all BLAST runs)...")
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
                                     list(foreground_file_id, background_file_id))
    } else {
      if (verbose) message("  - Extracting foreground GO terms via linkage (BLAST run ID: ", blast_param_id, ")...")
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
        WHERE c.file_id = ? AND r.file_id = ? AND br.blast_param_id = ?
      "
      foreground_go <- DBI::dbGetQuery(con, foreground_go_query,
                                     list(foreground_file_id, background_file_id, blast_param_id))
    }
  }

  # Background query with optional blast_param_id filtering
  if (is.null(blast_param_id)) {
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
    if (verbose) message("  - Extracting background GO terms (all BLAST runs)...")
    background_go <- DBI::dbGetQuery(con, background_go_query, list(background_file_id))
  } else {
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
      WHERE v.file_id = ? AND br.blast_param_id = ?
    "
    if (verbose) message("  - Extracting background GO terms (BLAST run ID: ", blast_param_id, ")...")
    background_go <- DBI::dbGetQuery(con, background_go_query, list(background_file_id, blast_param_id))
  }

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
#' @param significance_threshold Numeric. FDR threshold for significance classification. Default is 0.05
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
#' 
#' # Standard analysis (FDR < 0.05)
#' bp_results <- perform_go_enrichment(go_data, "BP")
#' 
#' # More lenient threshold (FDR < 0.1)
#' bp_results_lenient <- perform_go_enrichment(go_data, "BP", significance_threshold = 0.1)
#' 
#' head(bp_results)
#' }
#'
#' @keywords internal
perform_go_enrichment <- function(go_data, ontology = "BP", min_genes = 5, max_genes = 500, significance_threshold = 0.05, method = "clusterprofiler", verbose = TRUE) {

  if (verbose) message("Performing GO enrichment analysis for ontology: ", ontology, " using ", method, " method")

  # Route to appropriate enrichment method
  if (method == "clusterprofiler") {
    return(.perform_clusterprofiler_enrichment(go_data, ontology, min_genes, max_genes, significance_threshold, verbose))
  } else if (method == "legacy") {
    return(.perform_legacy_enrichment(go_data, ontology, min_genes, max_genes, significance_threshold, verbose))
  } else {
    stop("Invalid method. Must be 'clusterprofiler' or 'legacy'")
  }
}

#' Perform GO enrichment using clusterProfiler
#' @keywords internal
.perform_clusterprofiler_enrichment <- function(go_data, ontology, min_genes, max_genes, significance_threshold, verbose) {
  
  # Check if clusterProfiler is available
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("clusterProfiler package is required. Install with: BiocManager::install('clusterProfiler')")
  }
  
  # Convert go_data to clusterProfiler format
  clusterprofiler_data <- .convert_go_data_to_clusterprofiler(go_data, ontology, verbose)
  
  if (nrow(clusterprofiler_data$term2gene) == 0) {
    if (verbose) message("  - No ", ontology, " terms found for clusterProfiler analysis")
    return(data.frame())
  }
  
  # Run clusterProfiler enrichment
  tryCatch({
    enrichment_result <- clusterProfiler::enricher(
      gene = go_data$foreground$genes,
      universe = go_data$background$genes,
      TERM2GENE = clusterprofiler_data$term2gene,
      TERM2NAME = clusterprofiler_data$term2name,
      pvalueCutoff = 1.0,  # Get all results, filter later
      pAdjustMethod = "BH",
      minGSSize = min_genes,
      maxGSSize = max_genes
    )
    
    # Convert back to funseqR format
    return(.convert_clusterprofiler_to_funseqr(enrichment_result, ontology, significance_threshold, verbose))
    
  }, error = function(e) {
    warning("clusterProfiler enrichment failed: ", e$message)
    if (verbose) message("  - Falling back to legacy method")
    return(.perform_legacy_enrichment(go_data, ontology, min_genes, max_genes, significance_threshold, verbose))
  })
}

#' Perform GO enrichment using legacy hypergeometric method
#' @keywords internal
.perform_legacy_enrichment <- function(go_data, ontology, min_genes, max_genes, significance_threshold, verbose) {
  
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
    enrichment_results$p_adjusted < (significance_threshold / 5), "highly_significant",
    ifelse(enrichment_results$p_adjusted < significance_threshold, "significant",
           ifelse(enrichment_results$p_adjusted < (significance_threshold * 2), "trending", "not_significant"))
  )

  # Sort by significance and fold enrichment
  enrichment_results <- enrichment_results[order(enrichment_results$p_adjusted, -enrichment_results$fold_enrichment), ]

  if (verbose) {
    sig_count <- sum(enrichment_results$p_adjusted < significance_threshold)
    message("  - Enrichment analysis complete: ", nrow(enrichment_results), " terms tested, ", sig_count, " significantly enriched (FDR < ", significance_threshold, ")")
  }

  return(enrichment_results)
}

#' Convert go_data to clusterProfiler format
#' @keywords internal
.convert_go_data_to_clusterprofiler <- function(go_data, ontology, verbose) {
  
  # Map ontology codes
  ontology_map <- c("BP" = "P", "MF" = "F", "CC" = "C")
  category_code <- ontology_map[ontology]
  
  # Filter GO terms by category
  ontology_terms <- go_data$all_go_terms[go_data$all_go_terms$go_category == category_code, ]
  
  # Check for NULL or empty result
  if (is.null(ontology_terms) || nrow(ontology_terms) == 0) {
    return(list(term2gene = data.frame(term = character(0), gene = character(0)),
                term2name = data.frame(term = character(0), name = character(0))))
  }
  
  # Create TERM2GENE mapping (GO term -> gene)
  # For each gene, find which GO terms it has
  term2gene_list <- list()
  
  for (gene in go_data$background$genes) {
    if (gene %in% names(go_data$background$gene2go)) {
      gene_terms <- go_data$background$gene2go[[gene]]
      ontology_gene_terms <- intersect(gene_terms, ontology_terms$go_id)
      
      if (length(ontology_gene_terms) > 0) {
        for (term in ontology_gene_terms) {
          term2gene_list[[length(term2gene_list) + 1]] <- data.frame(
            term = term,
            gene = gene,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  if (length(term2gene_list) > 0) {
    term2gene <- do.call(rbind, term2gene_list)
  } else {
    term2gene <- data.frame(term = character(0), gene = character(0), stringsAsFactors = FALSE)
  }
  
  # Create TERM2NAME mapping (GO term -> term name)
  term2name <- unique(ontology_terms[, c("go_id", "go_term")])
  colnames(term2name) <- c("term", "name")
  
  if (verbose) {
    message("  - Converted ", nrow(term2name), " ", ontology, " terms for clusterProfiler")
    message("  - Total gene-term associations: ", nrow(term2gene))
  }
  
  return(list(
    term2gene = term2gene,
    term2name = term2name
  ))
}

#' Convert clusterProfiler results to funseqR format
#' @keywords internal
.convert_clusterprofiler_to_funseqr <- function(clusterprofiler_result, ontology, significance_threshold, verbose) {
  
  if (is.null(clusterprofiler_result) || nrow(clusterprofiler_result@result) == 0) {
    if (verbose) message("  - No terms found by clusterProfiler")
    return(data.frame())
  }
  
  cp_df <- clusterprofiler_result@result
  
  # Simple validation: check required clusterProfiler columns exist
  required_cols <- c("ID", "Description", "pvalue", "p.adjust")
  missing <- setdiff(required_cols, colnames(cp_df))
  if (length(missing) > 0) {
    stop("Missing clusterProfiler columns: ", paste(missing, collapse = ", "))
  }
  
  # Extract numeric values from ratios for calculations
  fg_count <- cp_df$Count
  bg_count <- as.numeric(sub("/.*", "", cp_df$BgRatio))
  total_fg <- as.numeric(sub(".*/", "", cp_df$GeneRatio))
  total_bg <- as.numeric(sub(".*/", "", cp_df$BgRatio))
  
  # Calculate expected count and fold enrichment correctly
  expected_count <- (bg_count / total_bg) * total_fg
  fold_enrichment <- ifelse(expected_count > 0, fg_count / expected_count, Inf)
  
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
    p_adjusted = as.numeric(cp_df$p.adjust),
    significance_level = ifelse(cp_df$p.adjust < (significance_threshold / 5), "highly_significant",
                               ifelse(cp_df$p.adjust < significance_threshold, "significant",
                                     ifelse(cp_df$p.adjust < (significance_threshold * 2), "trending", "not_significant"))),
    gene_ratio = cp_df$GeneRatio,
    bg_ratio = cp_df$BgRatio,
    qvalue = as.numeric(cp_df$qvalue),
    gene_ids = cp_df$geneID,
    stringsAsFactors = FALSE
  )
  
  if (verbose) {
    sig_count <- sum(funseqr_results$p_adjusted < significance_threshold, na.rm = TRUE)
    message("  - clusterProfiler analysis complete: ", nrow(funseqr_results), " terms tested, ", sig_count, " significantly enriched (FDR < ", significance_threshold, ")")
  }
  
  return(funseqr_results)
}

#' Extract KEGG pathways for foreground and background gene sets
#'
#' @param con Database connection object
#' @param foreground_file_id Integer. File ID of candidate/foreground dataset
#' @param background_file_id Integer. File ID of background dataset
#' @param blast_param_id Integer. Optional. Specific BLAST run ID to use for both datasets.
#'   If NULL, uses all available annotations. Default is NULL.
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return List containing foreground and background KEGG pathway data
#'
#' @keywords internal
extract_kegg_terms_for_enrichment <- function(con, foreground_file_id, background_file_id, 
                                            blast_param_id = NULL, verbose = TRUE) {

  if (verbose) message("Extracting KEGG pathways for enrichment analysis...")

  # Handle stored candidates case
  if (foreground_file_id == "stored") {
    if (verbose) message("  - Using stored candidate loci from database")
    
    # Check if candidates exist
    candidate_count <- DBI::dbGetQuery(con, "SELECT COUNT(*) as count FROM candidate_loci")$count
    if (candidate_count == 0) {
      stop("No candidate loci found in database. Run define_locus_statistics() or define_candidate_loci() first.")
    }
    
    # Extract KEGG annotations for stored candidates
    result <- .extract_stored_candidate_kegg_annotations(con, background_file_id, blast_param_id, verbose)
    return(result)
  }

  # First check if foreground file has direct annotations
  if (is.null(blast_param_id)) {
    foreground_direct_count <- DBI::dbGetQuery(con, "
      SELECT COUNT(*) as count FROM vcf_data v
      JOIN flanking_sequences fs ON v.vcf_id = fs.vcf_id
      JOIN blast_results br ON fs.flanking_id = br.flanking_id
      JOIN annotations a ON br.blast_result_id = a.blast_result_id
      WHERE v.file_id = ?
    ", list(foreground_file_id))$count
  } else {
    foreground_direct_count <- DBI::dbGetQuery(con, "
      SELECT COUNT(*) as count FROM vcf_data v
      JOIN flanking_sequences fs ON v.vcf_id = fs.vcf_id
      JOIN blast_results br ON fs.flanking_id = br.flanking_id
      JOIN annotations a ON br.blast_result_id = a.blast_result_id
      WHERE v.file_id = ? AND br.blast_param_id = ?
    ", list(foreground_file_id, blast_param_id))$count
  }

  if (foreground_direct_count > 0) {
    # Standard query for datasets with their own annotations
    if (is.null(blast_param_id)) {
      kegg_query <- "
        SELECT DISTINCT
          a.uniprot_accession,
          kr.kegg_id,
          kr.pathway_name
        FROM vcf_data v
        JOIN flanking_sequences fs ON v.vcf_id = fs.vcf_id
        JOIN blast_results br ON fs.flanking_id = br.flanking_id
        JOIN annotations a ON br.blast_result_id = a.blast_result_id
        JOIN kegg_references kr ON a.annotation_id = kr.annotation_id
        WHERE v.file_id = ?
      "
      if (verbose) message("  - Extracting foreground KEGG pathways (all BLAST runs)...")
      foreground_kegg <- DBI::dbGetQuery(con, kegg_query, list(foreground_file_id))
    } else {
      kegg_query <- "
        SELECT DISTINCT
          a.uniprot_accession,
          kr.kegg_id,
          kr.pathway_name
        FROM vcf_data v
        JOIN flanking_sequences fs ON v.vcf_id = fs.vcf_id
        JOIN blast_results br ON fs.flanking_id = br.flanking_id
        JOIN annotations a ON br.blast_result_id = a.blast_result_id
        JOIN kegg_references kr ON a.annotation_id = kr.annotation_id
        WHERE v.file_id = ? AND br.blast_param_id = ?
      "
      if (verbose) message("  - Extracting foreground KEGG pathways (BLAST run ID: ", blast_param_id, ")...")
      foreground_kegg <- DBI::dbGetQuery(con, kegg_query, list(foreground_file_id, blast_param_id))
    }
  } else {
    # For candidate files, use the linked annotations
    if (is.null(blast_param_id)) {
      if (verbose) message("  - Extracting foreground KEGG pathways via linkage (all BLAST runs)...")
      foreground_kegg_query <- "
        SELECT DISTINCT
          a.uniprot_accession,
          kr.kegg_id,
          kr.pathway_name
        FROM vcf_data c
        JOIN vcf_data r ON (c.chromosome = r.chromosome
                           AND c.position = r.position)
        JOIN flanking_sequences fs ON r.vcf_id = fs.vcf_id
        JOIN blast_results br ON fs.flanking_id = br.flanking_id
        JOIN annotations a ON br.blast_result_id = a.blast_result_id
        JOIN kegg_references kr ON a.annotation_id = kr.annotation_id
        WHERE c.file_id = ? AND r.file_id = ?
      "
      foreground_kegg <- DBI::dbGetQuery(con, foreground_kegg_query,
                                       list(foreground_file_id, background_file_id))
    } else {
      if (verbose) message("  - Extracting foreground KEGG pathways via linkage (BLAST run ID: ", blast_param_id, ")...")
      foreground_kegg_query <- "
        SELECT DISTINCT
          a.uniprot_accession,
          kr.kegg_id,
          kr.pathway_name
        FROM vcf_data c
        JOIN vcf_data r ON (c.chromosome = r.chromosome
                           AND c.position = r.position)
        JOIN flanking_sequences fs ON r.vcf_id = fs.vcf_id
        JOIN blast_results br ON fs.flanking_id = br.flanking_id
        JOIN annotations a ON br.blast_result_id = a.blast_result_id
        JOIN kegg_references kr ON a.annotation_id = kr.annotation_id
        WHERE c.file_id = ? AND r.file_id = ? AND br.blast_param_id = ?
      "
      foreground_kegg <- DBI::dbGetQuery(con, foreground_kegg_query,
                                       list(foreground_file_id, background_file_id, blast_param_id))
    }
  }

  # Background query with optional blast_param_id filtering
  if (is.null(blast_param_id)) {
    background_kegg_query <- "
      SELECT DISTINCT
        a.uniprot_accession,
        kr.kegg_id,
        kr.pathway_name
      FROM vcf_data v
      JOIN flanking_sequences fs ON v.vcf_id = fs.vcf_id
      JOIN blast_results br ON fs.flanking_id = br.flanking_id
      JOIN annotations a ON br.blast_result_id = a.blast_result_id
      JOIN kegg_references kr ON a.annotation_id = kr.annotation_id
      WHERE v.file_id = ?
    "
    if (verbose) message("  - Extracting background KEGG pathways (all BLAST runs)...")
    background_kegg <- DBI::dbGetQuery(con, background_kegg_query, list(background_file_id))
  } else {
    background_kegg_query <- "
      SELECT DISTINCT
        a.uniprot_accession,
        kr.kegg_id,
        kr.pathway_name
      FROM vcf_data v
      JOIN flanking_sequences fs ON v.vcf_id = fs.vcf_id
      JOIN blast_results br ON fs.flanking_id = br.flanking_id
      JOIN annotations a ON br.blast_result_id = a.blast_result_id
      JOIN kegg_references kr ON a.annotation_id = kr.annotation_id
      WHERE v.file_id = ? AND br.blast_param_id = ?
    "
    if (verbose) message("  - Extracting background KEGG pathways (BLAST run ID: ", blast_param_id, ")...")
    background_kegg <- DBI::dbGetQuery(con, background_kegg_query, list(background_file_id, blast_param_id))
  }

  # Create gene-to-pathway mapping lists
  foreground_gene2pathway <- split(foreground_kegg$kegg_id, foreground_kegg$uniprot_accession)
  background_gene2pathway <- split(background_kegg$kegg_id, background_kegg$uniprot_accession)

  # Get unique pathways and their details
  all_pathways <- rbind(
    foreground_kegg[, c("kegg_id", "pathway_name")],
    background_kegg[, c("kegg_id", "pathway_name")]
  )
  all_pathways <- all_pathways[!duplicated(all_pathways), ]

  if (verbose) {
    message("KEGG pathway extraction complete:")
    message("  - Foreground genes: ", length(unique(foreground_kegg$uniprot_accession)))
    message("  - Background genes: ", length(unique(background_kegg$uniprot_accession)))
    message("  - Total unique pathways: ", nrow(all_pathways))
  }

  return(list(
    foreground = list(
      genes = unique(foreground_kegg$uniprot_accession),
      gene2pathway = foreground_gene2pathway,
      pathways = foreground_kegg
    ),
    background = list(
      genes = unique(background_kegg$uniprot_accession),
      gene2pathway = background_gene2pathway,
      pathways = background_kegg
    ),
    all_pathways = all_pathways
  ))
}

#' Perform KEGG pathway enrichment analysis
#'
#' @param kegg_data List. Output from extract_kegg_terms_for_enrichment()
#' @param min_genes Integer. Minimum genes required for a pathway to be tested. Default is 5
#' @param max_genes Integer. Maximum genes for a pathway (to exclude very broad pathways). Default is 500
#' @param significance_threshold Numeric. FDR threshold for significance classification. Default is 0.05
#' @param method Character. Enrichment method: "clusterprofiler" or "legacy". Default is "clusterprofiler"
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return Data frame with enrichment results, sorted by adjusted p-value
#'
#' @keywords internal
perform_kegg_enrichment <- function(kegg_data, min_genes = 5, max_genes = 500, significance_threshold = 0.05, method = "clusterprofiler", verbose = TRUE) {

  if (verbose) message("Performing KEGG pathway enrichment analysis using ", method, " method")

  # Route to appropriate enrichment method
  if (method == "clusterprofiler") {
    return(.perform_clusterprofiler_kegg_enrichment(kegg_data, min_genes, max_genes, significance_threshold, verbose))
  } else if (method == "legacy") {
    return(.perform_legacy_kegg_enrichment(kegg_data, min_genes, max_genes, significance_threshold, verbose))
  } else {
    stop("Invalid method. Must be 'clusterprofiler' or 'legacy'")
  }
}

#' Perform KEGG enrichment using clusterProfiler
#' @keywords internal
.perform_clusterprofiler_kegg_enrichment <- function(kegg_data, min_genes, max_genes, significance_threshold, verbose) {
  
  # Check if clusterProfiler is available
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("clusterProfiler package is required. Install with: BiocManager::install('clusterProfiler')")
  }
  
  # Convert kegg_data to clusterProfiler format
  clusterprofiler_data <- .convert_kegg_data_to_clusterprofiler(kegg_data, verbose)
  
  if (nrow(clusterprofiler_data$term2gene) == 0) {
    if (verbose) message("  - No KEGG pathways found for clusterProfiler analysis")
    return(data.frame())
  }
  
  # Run clusterProfiler enrichment
  tryCatch({
    enrichment_result <- clusterProfiler::enricher(
      gene = kegg_data$foreground$genes,
      universe = kegg_data$background$genes,
      TERM2GENE = clusterprofiler_data$term2gene,
      TERM2NAME = clusterprofiler_data$term2name,
      pvalueCutoff = 1.0,  # Get all results, filter later
      pAdjustMethod = "BH",
      minGSSize = min_genes,
      maxGSSize = max_genes
    )
    
    # Convert back to funseqR format
    return(.convert_clusterprofiler_kegg_to_funseqr(enrichment_result, significance_threshold, verbose))
    
  }, error = function(e) {
    warning("clusterProfiler KEGG enrichment failed: ", e$message)
    if (verbose) message("  - Falling back to legacy method")
    return(.perform_legacy_kegg_enrichment(kegg_data, min_genes, max_genes, significance_threshold, verbose))
  })
}

#' Convert kegg_data to clusterProfiler format
#' @keywords internal
.convert_kegg_data_to_clusterprofiler <- function(kegg_data, verbose) {
  
  # Filter pathways
  pathways <- kegg_data$all_pathways
  
  # Check for NULL or empty pathways
  if (is.null(pathways) || nrow(pathways) == 0) {
    return(list(term2gene = data.frame(term = character(0), gene = character(0)),
                term2name = data.frame(term = character(0), name = character(0))))
  }
  
  # Create TERM2GENE mapping (pathway -> gene)
  term2gene_list <- list()
  
  for (gene in kegg_data$background$genes) {
    if (gene %in% names(kegg_data$background$gene2pathway)) {
      gene_pathways <- kegg_data$background$gene2pathway[[gene]]
      pathway_gene_pathways <- intersect(gene_pathways, pathways$kegg_id)
      
      if (length(pathway_gene_pathways) > 0) {
        for (pathway in pathway_gene_pathways) {
          term2gene_list[[length(term2gene_list) + 1]] <- data.frame(
            term = pathway,
            gene = gene,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  if (length(term2gene_list) > 0) {
    term2gene <- do.call(rbind, term2gene_list)
  } else {
    term2gene <- data.frame(term = character(0), gene = character(0), stringsAsFactors = FALSE)
  }
  
  # Create TERM2NAME mapping (pathway -> pathway name)
  term2name <- unique(pathways[, c("kegg_id", "pathway_name")])
  colnames(term2name) <- c("term", "name")
  
  if (verbose) {
    message("  - Converted ", nrow(term2name), " KEGG pathways for clusterProfiler")
    message("  - Total gene-pathway associations: ", nrow(term2gene))
  }
  
  return(list(
    term2gene = term2gene,
    term2name = term2name
  ))
}

#' Convert clusterProfiler KEGG results to funseqR format
#' @keywords internal
.convert_clusterprofiler_kegg_to_funseqr <- function(clusterprofiler_result, significance_threshold, verbose) {
  
  if (is.null(clusterprofiler_result) || nrow(clusterprofiler_result@result) == 0) {
    if (verbose) message("  - No pathways found by clusterProfiler")
    return(data.frame())
  }
  
  cp_df <- clusterprofiler_result@result
  
  # Convert to funseqR format with additional clusterProfiler fields
  funseqr_results <- data.frame(
    pathway_id = cp_df$ID,
    pathway_name = cp_df$Description,
    foreground_count = cp_df$Count,
    background_count = as.numeric(sub("/.*", "", cp_df$BgRatio)),
    total_foreground = as.numeric(sub(".*/", "", cp_df$GeneRatio)),
    total_background = as.numeric(sub(".*/", "", cp_df$BgRatio)),
    expected_count = cp_df$Count / cp_df$pvalue,  # Approximate
    fold_enrichment = cp_df$Count / (as.numeric(sub("/.*", "", cp_df$BgRatio)) / as.numeric(sub(".*/", "", cp_df$BgRatio)) * as.numeric(sub(".*/", "", cp_df$GeneRatio))),
    p_value = cp_df$pvalue,
    p_adjusted = cp_df$p.adjust,
    significance_level = ifelse(cp_df$p.adjust < (significance_threshold / 5), "highly_significant",
                               ifelse(cp_df$p.adjust < significance_threshold, "significant",
                                     ifelse(cp_df$p.adjust < (significance_threshold * 2), "trending", "not_significant"))),
    gene_ratio = cp_df$GeneRatio,
    bg_ratio = cp_df$BgRatio,
    qvalue = cp_df$qvalue,
    gene_ids = cp_df$geneID,
    stringsAsFactors = FALSE
  )
  
  if (verbose) {
    sig_count <- sum(funseqr_results$p_adjusted < significance_threshold, na.rm = TRUE)
    message("  - clusterProfiler KEGG analysis complete: ", nrow(funseqr_results), " pathways tested, ", sig_count, " significantly enriched (FDR < ", significance_threshold, ")")
  }
  
  return(funseqr_results)
}




# INTERNAL

#' Store ORA analysis results in the database
#'
#' @param con Database connection object
#' @param foreground_file_id Integer. File ID of foreground dataset
#' @param background_file_id Integer. File ID of background dataset
#' @param enrichment_results Data frame. Results from perform_go_enrichment() or perform_kegg_enrichment()
#' @param annotation_type Character. Type of annotation: "GO" or "KEGG"
#' @param term_type Character. Term type: "BP", "MF", "CC" for GO or "PATHWAY" for KEGG
#' @param parameters List. Analysis parameters for reproducibility
#' @param method Character. Enrichment method used. Default is "clusterprofiler"
#' @param blast_param_id Integer. Optional BLAST run ID used. Default is NULL
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return Integer. The analysis_id of the stored analysis
#'
store_ora_results <- function(con, foreground_file_id, background_file_id,
                             enrichment_results, annotation_type, term_type, parameters = NULL, method = "clusterprofiler", 
                             blast_param_id = NULL, verbose = TRUE) {

  if (verbose) message("Storing ", annotation_type, " ", term_type, " enrichment results in database...")

  # Ensure tables exist
  tables <- DBI::dbListTables(con)
  if (!all(c("ora_analyses", "ora_results") %in% tables)) {
    stop("ORA tables not found in database. Please upgrade schema.")
  }

  # Prepare analysis parameters
  if (is.null(parameters)) {
    parameters <- list()
  }
  params_json <- jsonlite::toJSON(parameters, auto_unbox = TRUE)

  # Insert analysis record
  analysis_query <- "
    INSERT INTO ora_analyses
    (foreground_file_id, background_file_id, blast_param_id, annotation_type, term_type, analysis_date,
     total_foreground_genes, total_background_genes, analysis_parameters, enrichment_method)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
  "

  DBI::dbExecute(con, analysis_query, list(
    foreground_file_id,
    background_file_id,
    blast_param_id,
    annotation_type,
    term_type,
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    if(nrow(enrichment_results) > 0) enrichment_results$total_foreground[1] else 0,
    if(nrow(enrichment_results) > 0) enrichment_results$total_background[1] else 0,
    params_json,
    method
  ))

  # Get the analysis ID
  analysis_id <- DBI::dbGetQuery(con, "SELECT last_insert_rowid() as id")$id

  # Insert results if any
  if (nrow(enrichment_results) > 0) {

    # Determine if this is clusterProfiler results (has additional columns)
    has_clusterprofiler_cols <- all(c("gene_ratio", "bg_ratio", "qvalue", "gene_ids") %in% colnames(enrichment_results))
    
    # Map column names based on annotation type
    if (annotation_type == "GO") {
      term_id_col <- "go_id"
      term_name_col <- "go_term"
    } else {
      term_id_col <- "pathway_id" 
      term_name_col <- "pathway_name"
    }
    
    # Debug: Check data types before storage
    if (verbose && nrow(enrichment_results) > 0) {
      sample_row <- enrichment_results[1, ]
      message("DEBUG: Storage data types - p_value: ", class(sample_row$p_value), 
              ", p_adjusted: ", class(sample_row$p_adjusted))
      message("DEBUG: Sample p_adjusted value: ", sample_row$p_adjusted)
    }
    
    # Simple row-by-row insertion with explicit column names
    for (i in 1:nrow(enrichment_results)) {
      row <- enrichment_results[i, ]
      
      if (has_clusterprofiler_cols) {
        # Insert with all clusterProfiler fields using explicit column names
        DBI::dbExecute(con, "
          INSERT INTO ora_results 
          (analysis_id, term_id, term_name, annotation_type, term_type,
           foreground_count, background_count, total_foreground, total_background,
           expected_count, fold_enrichment, p_value, p_adjusted, significance_level,
           gene_ratio, bg_ratio, qvalue, gene_ids)
          VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
          list(
            analysis_id,
            as.character(row[[term_id_col]]),
            as.character(row[[term_name_col]]),
            annotation_type,
            term_type,
            as.integer(row$foreground_count),
            as.integer(row$background_count),
            as.integer(row$total_foreground),
            as.integer(row$total_background),
            as.numeric(row$expected_count),
            as.numeric(row$fold_enrichment),
            as.numeric(row$p_value),
            as.numeric(row$p_adjusted),    # Explicit numeric conversion
            as.character(row$significance_level),
            as.character(row$gene_ratio),
            as.character(row$bg_ratio),
            as.numeric(row$qvalue),
            as.character(row$gene_ids)
          )
        )
      } else {
        # Legacy format without clusterProfiler extra fields
        DBI::dbExecute(con, "
          INSERT INTO ora_results 
          (analysis_id, term_id, term_name, annotation_type, term_type,
           foreground_count, background_count, total_foreground, total_background,
           expected_count, fold_enrichment, p_value, p_adjusted, significance_level)
          VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
          list(
            analysis_id,
            as.character(row[[term_id_col]]),
            as.character(row[[term_name_col]]),
            annotation_type,
            term_type,
            as.integer(row$foreground_count),
            as.integer(row$background_count),
            as.integer(row$total_foreground),
            as.integer(row$total_background),
            as.numeric(row$expected_count),
            as.numeric(row$fold_enrichment),
            as.numeric(row$p_value),
            as.numeric(row$p_adjusted),    # Explicit numeric conversion
            as.character(row$significance_level)
          )
        )
      }
    }
  }

  if (verbose) {
    message("Storage complete:")
    message("  - Analysis ID: ", analysis_id)
    message("  - Results stored: ", nrow(enrichment_results))
  }

  return(analysis_id)
}

#' Export GO terms for ReviGO summarization
#'
#' Export GO terms in ReviGO-compatible format for functional summarization using
#' the ReviGO web server (http://revigo.irb.hr). Supports three source types:
#' all annotations, candidate loci, and enriched terms.
#'
#' @param con Database connection object
#' @param source_type Character. Source of GO terms: "all_annotations", "candidate_loci", or "enriched_loci"
#' @param analysis_ids Integer vector. ORA analysis IDs to export (required for "enriched_loci")
#' @param candidate_loci Character, data.frame, or file path. Candidate loci specification (required for "candidate_loci")
#' @param vcf_file_id Integer. VCF file ID for automatic y-value extraction (optional for "candidate_loci")
#' @param full_y_values Numeric vector. Full y-values for all SNPs in VCF order (optional for "candidate_loci")
#' @param candidate_locus_pvalues Numeric vector. P-values for each candidate locus, same order as candidate_loci (optional for "candidate_loci")
#' @param include_pvalues Logical. Whether to include p-values in output (for "candidate_loci"). Default is TRUE
#' @param significance_threshold Numeric. P-value threshold for enriched terms (for "enriched_loci"). Default is 0.05
#' @param ontologies Character vector. GO ontologies to include: c("BP", "MF", "CC"). If NULL, includes all
#' @param output_file Character. Output file path. If NULL, returns data.frame instead of writing file
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return If output_file is NULL, returns data.frame. Otherwise returns file path and writes ReviGO-compatible file
#'
#' @details
#' This function exports GO terms for summarization using the widely-used ReviGO web server.
#' Three source types are supported:
#'
#' \\strong{All Annotations} (\\code{source_type = "all_annotations"}):
#' Exports all GO terms found in the database. Output contains only GO term IDs (no p-values).
#' Useful for summarizing the complete functional landscape of your dataset.
#'
#' \\strong{Candidate Loci} (\\code{source_type = "candidate_loci"}):
#' Exports GO terms associated with user-specified candidate loci. Two output modes:
#' \\itemize{
#'   \\item \\code{include_pvalues = FALSE}: GO term IDs only
#'   \\item \\code{include_pvalues = TRUE}: GO terms with p-values derived from statistical evidence
#' }
#' When p-values are included, each GO term receives the minimum p-value among all associated loci,
#' effectively weighting functional terms by the strength of statistical evidence for selection.
#' 
#' P-values can be provided in two ways:
#' \\itemize{
#'   \\item \\strong{Automatic extraction}: Provide \\code{vcf_file_id} and \\code{full_y_values} to automatically extract p-values for candidate loci from the full dataset
#'   \\item \\strong{Direct specification}: Provide \\code{candidate_locus_pvalues} with p-values matching candidate loci order
#' }
#'
#' \\strong{Enriched Loci} (\\code{source_type = "enriched_loci"}):
#' Exports statistically enriched GO terms from ORA analysis with their real statistical p-values.
#' Uses results from \\code{ora()} analysis.
#'
#' \\strong{Input Formats:}
#' The \\code{candidate_loci} parameter supports multiple input formats:
#' \\itemize{
#'   \\item VCF file path: "candidates.vcf"
#'   \\item VCF file ID: 1 (integer referencing database)
#'   \\item Coordinate strings: c("LG4:3814415", "LG8:22215908")
#'   \\item Data.frame: data.frame(chromosome = c("LG4", "LG8"), position = c(3814415, 22215908))
#'   \\item BED file path: "regions.bed"
#'   \\item \\strong{NULL}: Automatically uses stored candidate loci from \\code{define_locus_statistics()} or \\code{define_candidate_loci()}
#' }
#' 
#' \\strong{Streamlined Workflow Integration:}
#' This function integrates with the streamlined workflow:
#' \\itemize{
#'   \\item When \\code{candidate_loci = NULL}, automatically uses stored candidate loci from database
#'   \\item When \\code{include_pvalues = TRUE}, automatically uses stored statistics for p-values (if available)
#'   \\item Falls back to manual specification if stored data not available
#' }
#'
#' \\strong{Output Formats:}
#' ReviGO-compatible tab-separated format:
#' \\itemize{
#'   \\item With p-values: "\\% GOterm\\tenrichment_P-value" header + "GO:0022402\\t1.74E-15"
#'   \\item Without p-values: "GO:0022402" (one GO term per line)
#' }
#'
#' @examples
#' \\dontrun{
#' con <- connect_funseq_db("analysis.db")
#'
#' # Export all GO annotations (functional landscape summary)
#' export_revigo_input(con, source_type = "all_annotations",
#'                     output_file = "all_functions.txt")
#'
#' # Export candidate loci - automatic y-value extraction (recommended)
#' export_revigo_input(con, source_type = "candidate_loci",
#'                     candidate_loci = "candidates.vcf",
#'                     vcf_file_id = 1,                    # full dataset
#'                     full_y_values = my_rda_pvalues,    # all SNP p-values
#'                     include_pvalues = TRUE,
#'                     output_file = "candidate_functions.txt")
#'
#' # Export candidate loci - direct p-value specification
#' export_revigo_input(con, source_type = "candidate_loci",
#'                     candidate_loci = c("LG4:3814415", "LG8:22215908"),
#'                     candidate_locus_pvalues = c(0.001, 0.003),
#'                     include_pvalues = TRUE,
#'                     output_file = "candidate_functions.txt")
#'
#' # Export candidate loci functions without p-values
#' export_revigo_input(con, source_type = "candidate_loci",
#'                     candidate_loci = "candidates.vcf",
#'                     include_pvalues = FALSE,
#'                     output_file = "candidate_go_terms.txt")
#'
#' # Export enriched terms from ORA analysis
#' export_revigo_input(con, source_type = "enriched_loci",
#'                     analysis_ids = c(1, 2, 3),
#'                     significance_threshold = 0.01,
#'                     ontologies = c("BP", "MF"),
#'                     output_file = "enriched_functions.txt")
#' 
#' # Streamlined workflow with stored data
#' # Step 1: Define statistics and candidates
#' define_locus_statistics(con, my_pvalue_data, candidate_threshold = 0.01)
#' 
#' # Step 2: Export candidates to ReviGO (uses stored data automatically)
#' export_revigo_input(con, source_type = "candidate_loci",
#'                     output_file = "candidates_for_revigo.txt")
#' }
#'
#' @export
export_revigo_input <- function(con, source_type = c("all_annotations", "candidate_loci", "enriched_loci"),
                               analysis_ids = NULL, candidate_loci = NULL, 
                               vcf_file_id = NULL, full_y_values = NULL, candidate_locus_pvalues = NULL,
                               include_pvalues = TRUE, significance_threshold = 0.05,
                               ontologies = NULL, output_file = NULL, verbose = TRUE) {
  
  # Validate parameters
  source_type <- match.arg(source_type)
  
  if (verbose) message("=== Exporting GO terms for ReviGO ===")
  if (verbose) message("Source type: ", source_type)
  
  # Validate required parameters for each source type
  if (source_type == "enriched_loci" && is.null(analysis_ids)) {
    stop("analysis_ids is required for source_type = 'enriched_loci'")
  }
  
  if (source_type == "candidate_loci" && is.null(candidate_loci)) {
    # Check if stored candidate loci exist
    stored_candidates <- DBI::dbGetQuery(con, "SELECT COUNT(*) as count FROM candidate_loci")$count
    if (stored_candidates == 0) {
      stop("candidate_loci is required for source_type = 'candidate_loci' when no stored candidates exist. Use define_locus_statistics() or define_candidate_loci() first, or provide candidate_loci parameter.")
    } else {
      if (verbose) message("Using ", stored_candidates, " stored candidate loci from database")
      # Create candidate_loci specification from stored data
      candidate_loci <- "stored"  # Flag to use stored data
    }
  }
  
  # Validate p-value specification for candidate_loci
  if (source_type == "candidate_loci" && include_pvalues) {
    has_direct_pvalues <- !is.null(candidate_locus_pvalues)
    has_auto_extraction <- !is.null(vcf_file_id) && !is.null(full_y_values)
    
    if (!has_direct_pvalues && !has_auto_extraction) {
      stop("When source_type = 'candidate_loci' and include_pvalues = TRUE, you must provide either:\n",
           "  - candidate_locus_pvalues (direct specification), or\n",
           "  - vcf_file_id and full_y_values (automatic extraction)")
    }
    
    if (has_direct_pvalues && has_auto_extraction) {
      if (verbose) message("Both candidate_locus_pvalues and automatic extraction provided. Using candidate_locus_pvalues directly.")
    }
  }
  
  # Extract GO terms based on source type
  if (source_type == "all_annotations") {
    go_data <- .extract_all_go_annotations(con, ontologies, verbose)
    
  } else if (source_type == "candidate_loci") {
    go_data <- .extract_candidate_go_annotations(con, candidate_loci, vcf_file_id, full_y_values, 
                                                candidate_locus_pvalues, include_pvalues, ontologies, verbose)
    
  } else if (source_type == "enriched_loci") {
    go_data <- .extract_enriched_go_terms(con, analysis_ids, significance_threshold, 
                                         ontologies, verbose)
  }
  
  if (nrow(go_data) == 0) {
    if (verbose) message("No GO terms found matching criteria")
    if (!is.null(output_file)) {
      # Create empty file
      writeLines("", output_file)
      if (verbose) message("Empty file written to: ", output_file)
      return(output_file)
    } else {
      return(data.frame())
    }
  }
  
  # Format for ReviGO
  revigo_data <- .format_revigo_output(go_data, include_pvalues = "p_value" %in% names(go_data), verbose)
  
  # Output results
  if (!is.null(output_file)) {
    .write_revigo_file(revigo_data, output_file, include_pvalues = "p_value" %in% names(go_data), verbose)
    if (verbose) message("ReviGO file written to: ", output_file)
    return(output_file)
  } else {
    if (verbose) message("Returning ", nrow(revigo_data), " GO terms as data.frame")
    return(revigo_data)
  }
}

# Helper functions for export_revigo_input()

#' Extract all GO annotations from database
#' @keywords internal
.extract_all_go_annotations <- function(con, ontologies, verbose) {
  
  if (verbose) message("  - Extracting all GO annotations from database...")
  
  # Build ontology filter
  ontology_filter <- ""
  if (!is.null(ontologies)) {
    ontology_list <- paste0("'", ontologies, "'", collapse = ", ")
    ontology_filter <- paste0(" AND gt.ontology IN (", ontology_list, ")")
    if (verbose) message("    - Filtering to ontologies: ", paste(ontologies, collapse = ", "))
  }
  
  query <- paste0("
    SELECT DISTINCT gt.go_id as term_id
    FROM go_terms gt
    JOIN annotations a ON gt.annotation_id = a.annotation_id
    WHERE gt.go_id IS NOT NULL AND gt.go_id != ''",
    ontology_filter, "
    ORDER BY gt.go_id
  ")
  
  result <- DBI::dbGetQuery(con, query)
  
  if (verbose) message("    - Found ", nrow(result), " unique GO terms")
  
  return(result)
}

#' Extract GO annotations for candidate loci
#' @keywords internal
.extract_candidate_go_annotations <- function(con, candidate_loci, vcf_file_id, full_y_values,
                                            candidate_locus_pvalues, include_pvalues, ontologies, verbose) {
  
  if (verbose) message("  - Extracting GO annotations for candidate loci...")
  
  # Parse candidate loci using existing infrastructure
  if (is.character(candidate_loci) && length(candidate_loci) == 1) {
    # Check for stored candidates flag
    if (candidate_loci == "stored") {
      if (verbose) message("    - Using stored candidate loci from database")
      loci_coords <- DBI::dbGetQuery(con, "
        SELECT chromosome, position FROM candidate_loci ORDER BY candidate_id
      ")
    } else if (file.exists(candidate_loci)) {
      if (verbose) message("    - Processing file: ", candidate_loci)
      if (grepl("\\.bed$", candidate_loci, ignore.case = TRUE)) {
        # BED file
        loci_coords <- .parse_bed_file(candidate_loci, verbose)
      } else {
        # Assume VCF file - import and get coordinates
        temp_import <- import_vcf(con, candidate_loci)
        vcf_coords <- DBI::dbGetQuery(con, "
          SELECT chromosome, position FROM vcf_data 
          WHERE file_id = ? ORDER BY vcf_id
        ", list(temp_import$file_id))
        loci_coords <- vcf_coords
      }
    } else {
      # Try as VCF file ID
      vcf_id <- suppressWarnings(as.integer(candidate_loci))
      if (!is.na(vcf_id)) {
        if (verbose) message("    - Using VCF file ID: ", vcf_id)
        vcf_coords <- DBI::dbGetQuery(con, "
          SELECT chromosome, position FROM vcf_data 
          WHERE file_id = ? ORDER BY vcf_id
        ", list(vcf_id))
        loci_coords <- vcf_coords
      } else {
        # Single coordinate string
        loci_coords <- .parse_loci_input(candidate_loci, verbose)
      }
    }
  } else {
    # Multiple coordinates or data.frame
    loci_coords <- .parse_loci_input(candidate_loci, verbose)
  }
  
  if (nrow(loci_coords) == 0) {
    if (verbose) message("    - No valid candidate loci found")
    return(data.frame(term_id = character(0), p_value = numeric(0)))
  }
  
  # Extract y-values automatically if requested and not directly provided
  if (include_pvalues && is.null(candidate_locus_pvalues)) {
    # Try stored statistics first
    stored_stats <- DBI::dbGetQuery(con, "SELECT COUNT(*) as count FROM locus_statistics")$count
    if (stored_stats > 0) {
      if (verbose) message("    - Extracting y-values for candidate loci from stored statistics...")
      candidate_locus_pvalues <- .extract_candidate_yvalues_from_stored(con, loci_coords, verbose)
    } else if (!is.null(vcf_file_id) && !is.null(full_y_values)) {
      if (verbose) message("    - Extracting y-values for candidate loci from full dataset...")
      candidate_locus_pvalues <- .extract_candidate_yvalues(con, loci_coords, vcf_file_id, full_y_values, verbose)
    }
  }
  
  # Validate p-values length if provided
  if (include_pvalues && !is.null(candidate_locus_pvalues)) {
    if (length(candidate_locus_pvalues) != nrow(loci_coords)) {
      stop("Length of candidate_locus_pvalues (", length(candidate_locus_pvalues), 
           ") does not match number of loci (", nrow(loci_coords), ")")
    }
  }
  
  # Build ontology filter
  ontology_filter <- ""
  if (!is.null(ontologies)) {
    ontology_list <- paste0("'", ontologies, "'", collapse = ", ")
    ontology_filter <- paste0(" AND gt.ontology IN (", ontology_list, ")")
  }
  
  # Query GO terms for these loci
  go_results <- data.frame(term_id = character(0), locus_chromosome = character(0), 
                          locus_position = numeric(0), stringsAsFactors = FALSE)
  
  for (i in 1:nrow(loci_coords)) {
    chr <- loci_coords$chromosome[i]
    pos <- loci_coords$position[i]
    
    query <- paste0("
      SELECT DISTINCT 
        gt.go_id as term_id,
        vd.chromosome as locus_chromosome,
        vd.position as locus_position
      FROM vcf_data vd
      JOIN flanking_sequences fs ON vd.vcf_id = fs.vcf_id
      JOIN blast_results br ON fs.flanking_id = br.flanking_id
      JOIN annotations a ON br.blast_result_id = a.blast_result_id
      JOIN go_terms gt ON a.annotation_id = gt.annotation_id
      WHERE vd.chromosome = ? AND vd.position = ?
        AND gt.go_id IS NOT NULL AND gt.go_id != ''",
      ontology_filter
    )
    
    locus_results <- DBI::dbGetQuery(con, query, list(chr, pos))
    go_results <- rbind(go_results, locus_results)
  }
  
  if (nrow(go_results) == 0) {
    if (verbose) message("    - No GO annotations found for candidate loci")
    if (include_pvalues) {
      return(data.frame(term_id = character(0), p_value = numeric(0)))
    } else {
      return(data.frame(term_id = character(0)))
    }
  }
  
  # Aggregate by GO term and assign p-values if requested
  if (include_pvalues && !is.null(candidate_locus_pvalues)) {
    # Create lookup for locus p-values
    locus_pvalues <- data.frame(
      chromosome = loci_coords$chromosome,
      position = loci_coords$position,
      p_value = candidate_locus_pvalues,
      stringsAsFactors = FALSE
    )
    
    # Merge p-values with GO results
    go_results <- merge(go_results, locus_pvalues, 
                       by.x = c("locus_chromosome", "locus_position"),
                       by.y = c("chromosome", "position"))
    
    # Aggregate: take minimum p-value per GO term
    result <- aggregate(p_value ~ term_id, data = go_results, FUN = min)
    result <- result[order(result$p_value), ]
    
    if (verbose) message("    - Found ", nrow(result), " GO terms with p-values")
    
  } else {
    # Just unique GO terms
    result <- data.frame(term_id = unique(go_results$term_id), stringsAsFactors = FALSE)
    result <- result[order(result$term_id), , drop = FALSE]
    
    if (verbose) message("    - Found ", nrow(result), " unique GO terms")
  }
  
  return(result)
}

#' Extract enriched GO terms from ORA results
#' @keywords internal
.extract_enriched_go_terms <- function(con, analysis_ids, significance_threshold, ontologies, verbose) {
  
  if (verbose) message("  - Extracting enriched GO terms from ORA analysis...")
  
  # Use existing function from process_annotations.R
  enrichment_data <- .extract_enrichment_data(con, analysis_ids, significance_threshold, verbose = FALSE)
  
  if (nrow(enrichment_data) == 0) {
    if (verbose) message("    - No enriched terms found")
    return(data.frame(term_id = character(0), p_value = numeric(0)))
  }
  
  # Filter to GO terms only
  go_enrichment <- enrichment_data[enrichment_data$annotation_type == "GO", ]
  
  if (nrow(go_enrichment) == 0) {
    if (verbose) message("    - No enriched GO terms found")
    return(data.frame(term_id = character(0), p_value = numeric(0)))
  }
  
  # Filter by ontologies if specified
  if (!is.null(ontologies)) {
    go_enrichment <- go_enrichment[go_enrichment$term_type %in% ontologies, ]
    if (verbose) message("    - Filtering to ontologies: ", paste(ontologies, collapse = ", "))
  }
  
  if (nrow(go_enrichment) == 0) {
    if (verbose) message("    - No GO terms found after ontology filtering")
    return(data.frame(term_id = character(0), p_value = numeric(0)))
  }
  
  # Format results
  result <- data.frame(
    term_id = go_enrichment$term_id,
    p_value = as.numeric(go_enrichment$p_value),
    stringsAsFactors = FALSE
  )
  
  # Remove rows with missing p-values and sort
  result <- result[!is.na(result$p_value), ]
  result <- result[order(result$p_value), ]
  
  if (verbose) message("    - Found ", nrow(result), " enriched GO terms")
  
  return(result)
}

#' Format GO data for ReviGO output
#' @keywords internal
.format_revigo_output <- function(go_data, include_pvalues, verbose) {
  
  if (verbose) message("  - Formatting data for ReviGO...")
  
  if (include_pvalues && "p_value" %in% names(go_data)) {
    # Format with p-values in scientific notation
    result <- data.frame(
      term_id = go_data$term_id,
      p_value = sprintf("%.2E", go_data$p_value),
      stringsAsFactors = FALSE
    )
    if (verbose) message("    - Formatted ", nrow(result), " GO terms with p-values")
  } else {
    # Just GO terms
    result <- data.frame(
      term_id = go_data$term_id,
      stringsAsFactors = FALSE
    )
    if (verbose) message("    - Formatted ", nrow(result), " GO terms")
  }
  
  return(result)
}

#' Write ReviGO-compatible file
#' @keywords internal
.write_revigo_file <- function(revigo_data, output_file, include_pvalues, verbose) {
  
  if (verbose) message("  - Writing ReviGO file...")
  
  if (include_pvalues && ncol(revigo_data) > 1) {
    # Write with header and p-values
    lines <- c("% GOterm\tenrichment_P-value",
               paste(revigo_data$term_id, revigo_data$p_value, sep = "\t"))
  } else {
    # Write just GO terms
    lines <- revigo_data$term_id
  }
  
  writeLines(lines, output_file)
  
  if (verbose) message("    - Written ", length(lines) - ifelse(include_pvalues, 1, 0), " GO terms")
}

#' Extract y-values for candidate loci from full VCF dataset
#' @keywords internal
.extract_candidate_yvalues <- function(con, candidate_coords, vcf_file_id, full_y_values, verbose) {
  
  # Get full VCF coordinates in order (same approach as Manhattan plot)
  vcf_coords <- DBI::dbGetQuery(con, "
    SELECT vcf_id, chromosome, position, ref, alt
    FROM vcf_data
    WHERE file_id = ?
    ORDER BY vcf_id
  ", list(vcf_file_id))
  
  if (nrow(vcf_coords) == 0) {
    stop("No VCF data found for file_id: ", vcf_file_id)
  }
  
  # Validate y_values length (same validation as Manhattan plot)
  if (length(full_y_values) != nrow(vcf_coords)) {
    stop("Length mismatch: full_y_values has ", length(full_y_values),
         " values but VCF has ", nrow(vcf_coords), " variants")
  }
  
  if (verbose) message("      - Matching ", nrow(candidate_coords), " candidate loci against ", nrow(vcf_coords), " VCF variants")
  
  # Extract y-values for candidate loci
  candidate_pvalues <- numeric(nrow(candidate_coords))
  matched_count <- 0
  
  for (i in 1:nrow(candidate_coords)) {
    chr <- candidate_coords$chromosome[i]
    pos <- candidate_coords$position[i]
    
    # Find matching VCF entry
    match_idx <- which(vcf_coords$chromosome == chr & vcf_coords$position == pos)
    
    if (length(match_idx) > 0) {
      candidate_pvalues[i] <- full_y_values[match_idx[1]]
      matched_count <- matched_count + 1
    } else {
      if (verbose) message("      - Warning: Candidate locus not found in VCF: ", chr, ":", pos)
      candidate_pvalues[i] <- NA
    }
  }
  
  if (verbose) message("      - Successfully matched ", matched_count, "/", nrow(candidate_coords), " candidate loci")
  
  # Report NA values but don't remove them here - let the calling function handle it
  if (any(is.na(candidate_pvalues))) {
    na_count <- sum(is.na(candidate_pvalues))
    if (verbose) message("      - Warning: ", na_count, " candidate loci not found in VCF and will be excluded from p-value analysis")
  }
  
  return(candidate_pvalues)
}

# Helper function to extract y-values from stored statistics
.extract_candidate_yvalues_from_stored <- function(con, candidate_coords, verbose) {
  
  if (verbose) message("      - Matching candidate loci to stored statistics...")
  
  candidate_pvalues <- numeric(nrow(candidate_coords))
  matched_count <- 0
  
  for (i in 1:nrow(candidate_coords)) {
    chr <- candidate_coords$chromosome[i]
    pos <- candidate_coords$position[i]
    
    # Query stored statistics for this coordinate
    stat_result <- DBI::dbGetQuery(con, "
      SELECT statistic FROM locus_statistics 
      WHERE chromosome = ? AND position = ?
    ", list(chr, pos))
    
    if (nrow(stat_result) > 0) {
      candidate_pvalues[i] <- stat_result$statistic[1]
      matched_count <- matched_count + 1
    } else {
      if (verbose) message("      - Warning: Candidate locus not found in stored statistics: ", chr, ":", pos)
      candidate_pvalues[i] <- NA
    }
  }
  
  if (verbose) message("      - Successfully matched ", matched_count, "/", nrow(candidate_coords), " candidate loci")
  
  # Report NA values
  if (any(is.na(candidate_pvalues))) {
    na_count <- sum(is.na(candidate_pvalues))
    if (verbose) message("      - Warning: ", na_count, " candidate loci not found in stored statistics and will be excluded from p-value analysis")
  }
  
  return(candidate_pvalues)
}

#' Extract GO annotations for stored candidates
#' @keywords internal
.extract_stored_candidate_go_annotations <- function(con, background_file_id, blast_param_id = NULL, verbose = TRUE) {
  
  if (verbose) message("  - Extracting GO annotations for stored candidates...")
  
  # Build base conditions for BLAST filtering
  blast_condition <- if (!is.null(blast_param_id)) {
    "AND br.blast_param_id = ?"
  } else {
    ""
  }
  
  # Query for candidate (foreground) data - use stored candidates
  candidate_query <- paste0("
    SELECT DISTINCT 
      gt.go_id as term_id,
      a.uniprot_accession as gene_id,
      gt.go_term as term_name,
      gt.go_category as ontology,
      vd.chromosome,
      vd.position
    FROM candidate_loci cl
    JOIN vcf_data vd ON cl.chromosome = vd.chromosome AND cl.position = vd.position
    JOIN flanking_sequences fs ON vd.vcf_id = fs.vcf_id
    JOIN blast_results br ON fs.flanking_id = br.flanking_id ", blast_condition, "
    JOIN annotations a ON br.blast_result_id = a.blast_result_id
    JOIN go_terms gt ON a.annotation_id = gt.annotation_id
    WHERE vd.file_id = ?
    ORDER BY gt.go_id, a.uniprot_accession
  ")
  
  # Execute candidate query
  candidate_params <- if (!is.null(blast_param_id)) {
    list(blast_param_id, background_file_id)
  } else {
    list(background_file_id)
  }
  
  candidate_data <- DBI::dbGetQuery(con, candidate_query, candidate_params)
  
  # Query for background data (same as original function)
  background_query <- paste0("
    SELECT DISTINCT 
      gt.go_id as term_id,
      a.uniprot_accession as gene_id,
      gt.go_term as term_name,
      gt.go_category as ontology,
      vd.chromosome,
      vd.position
    FROM vcf_data vd
    JOIN flanking_sequences fs ON vd.vcf_id = fs.vcf_id
    JOIN blast_results br ON fs.flanking_id = br.flanking_id ", blast_condition, "
    JOIN annotations a ON br.blast_result_id = a.blast_result_id
    JOIN go_terms gt ON a.annotation_id = gt.annotation_id
    WHERE vd.file_id = ?
    ORDER BY gt.go_id, a.uniprot_accession
  ")
  
  background_data <- DBI::dbGetQuery(con, background_query, candidate_params)
  
  if (verbose) {
    message("    - Found ", nrow(candidate_data), " candidate GO associations")
    message("    - Found ", nrow(background_data), " background GO associations")
  }
  
  # Format results same as original function
  if (nrow(candidate_data) == 0) {
    warning("No GO annotations found for stored candidate loci")
    return(list(
      foreground = list(genes = character(0), terms = character(0)),
      background = list(genes = character(0), terms = character(0)),
      term2gene_fg = data.frame(term = character(0), gene = character(0)),
      term2gene_bg = data.frame(term = character(0), gene = character(0)),
      gene2name = data.frame(gene = character(0), name = character(0))
    ))
  }
  
  # Process the data into the expected format
  fg_genes <- unique(candidate_data$gene_id)
  fg_terms <- unique(candidate_data$term_id)
  bg_genes <- unique(background_data$gene_id)
  bg_terms <- unique(background_data$term_id)
  
  term2gene_fg <- candidate_data[, c("term_id", "gene_id")]
  names(term2gene_fg) <- c("term", "gene")
  
  term2gene_bg <- background_data[, c("term_id", "gene_id")]
  names(term2gene_bg) <- c("term", "gene")
  
  # Create gene2name mapping (simplified)
  all_genes <- unique(c(candidate_data$gene_id, background_data$gene_id))
  gene2name <- data.frame(
    gene = all_genes,
    name = all_genes,  # Use gene ID as name for simplicity
    stringsAsFactors = FALSE
  )
  
  return(list(
    foreground = list(genes = fg_genes, terms = fg_terms),
    background = list(genes = bg_genes, terms = bg_terms),
    term2gene_fg = term2gene_fg,
    term2gene_bg = term2gene_bg,
    gene2name = gene2name
  ))
}

#' Extract KEGG annotations for stored candidates
#' @keywords internal
.extract_stored_candidate_kegg_annotations <- function(con, background_file_id, blast_param_id = NULL, verbose = TRUE) {
  
  if (verbose) message("  - Extracting KEGG annotations for stored candidates...")
  
  # Build base conditions for BLAST filtering
  blast_condition <- if (!is.null(blast_param_id)) {
    "AND br.blast_param_id = ?"
  } else {
    ""
  }
  
  # Query for candidate (foreground) data - use stored candidates
  candidate_query <- paste0("
    SELECT DISTINCT 
      kr.kegg_id as term_id,
      a.uniprot_accession as gene_id,
      kr.pathway_name as term_name,
      vd.chromosome,
      vd.position
    FROM candidate_loci cl
    JOIN vcf_data vd ON cl.chromosome = vd.chromosome AND cl.position = vd.position
    JOIN flanking_sequences fs ON vd.vcf_id = fs.vcf_id
    JOIN blast_results br ON fs.flanking_id = br.flanking_id ", blast_condition, "
    JOIN annotations a ON br.blast_result_id = a.blast_result_id
    JOIN kegg_references kr ON a.annotation_id = kr.annotation_id
    WHERE vd.file_id = ?
    ORDER BY kr.kegg_id, a.uniprot_accession
  ")
  
  # Execute candidate query
  candidate_params <- if (!is.null(blast_param_id)) {
    list(blast_param_id, background_file_id)
  } else {
    list(background_file_id)
  }
  
  candidate_data <- DBI::dbGetQuery(con, candidate_query, candidate_params)
  
  # Query for background data
  background_query <- paste0("
    SELECT DISTINCT 
      kr.kegg_id as term_id,
      a.uniprot_accession as gene_id,
      kr.pathway_name as term_name,
      vd.chromosome,
      vd.position
    FROM vcf_data vd
    JOIN flanking_sequences fs ON vd.vcf_id = fs.vcf_id
    JOIN blast_results br ON fs.flanking_id = br.flanking_id ", blast_condition, "
    JOIN annotations a ON br.blast_result_id = a.blast_result_id
    JOIN kegg_references kr ON a.annotation_id = kr.annotation_id
    WHERE vd.file_id = ?
    ORDER BY kr.kegg_id, a.uniprot_accession
  ")
  
  background_data <- DBI::dbGetQuery(con, background_query, candidate_params)
  
  if (verbose) {
    message("    - Found ", nrow(candidate_data), " candidate KEGG associations")
    message("    - Found ", nrow(background_data), " background KEGG associations")
  }
  
  # Format results same as original function
  if (nrow(candidate_data) == 0) {
    warning("No KEGG annotations found for stored candidate loci")
    return(list(
      foreground = list(genes = character(0), terms = character(0)),
      background = list(genes = character(0), terms = character(0)),
      term2gene_fg = data.frame(term = character(0), gene = character(0)),
      term2gene_bg = data.frame(term = character(0), gene = character(0)),
      gene2name = data.frame(gene = character(0), name = character(0))
    ))
  }
  
  # Process the data into the expected format
  fg_genes <- unique(candidate_data$gene_id)
  fg_terms <- unique(candidate_data$term_id)
  bg_genes <- unique(background_data$gene_id)
  bg_terms <- unique(background_data$term_id)
  
  term2gene_fg <- candidate_data[, c("term_id", "gene_id")]
  names(term2gene_fg) <- c("term", "gene")
  
  term2gene_bg <- background_data[, c("term_id", "gene_id")]
  names(term2gene_bg) <- c("term", "gene")
  
  # Create gene2name mapping (simplified)
  all_genes <- unique(c(candidate_data$gene_id, background_data$gene_id))
  gene2name <- data.frame(
    gene = all_genes,
    name = all_genes,  # Use gene ID as name for simplicity
    stringsAsFactors = FALSE
  )
  
  return(list(
    foreground = list(genes = fg_genes, terms = fg_terms),
    background = list(genes = bg_genes, terms = bg_terms),
    term2gene_fg = term2gene_fg,
    term2gene_bg = term2gene_bg,
    gene2name = gene2name
  ))
}

