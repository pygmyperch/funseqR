# EXPORTED

#' GO Enrichment Analysis Functions
#'
#' This file contains functions for performing Gene Ontology (GO) enrichment analysis
#' on candidate loci vs background datasets, with visualization capabilities.
#'

#' Import candidate adaptive loci and link to existing annotations
#'
#' @param con Database connection object
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
#' candidate_import <- import_candidate_loci(con, "candidates.vcf", 1)
#' }
#'
#' @export
import_candidate_loci <- function(con, candidate_vcf_file, background_file_id, verbose = TRUE) {

  if (verbose) message("Importing candidate loci from: ", candidate_vcf_file)

  # Import candidate VCF using existing function
  candidate_import <- import_vcf_to_db(con, candidate_vcf_file)

  if (verbose) message("Creating BED file for candidate loci...")

  # Create BED file for candidate loci
  candidate_bed <- vcf2bed_db(con, candidate_import$file_id)

  if (verbose) message("Linking candidates to existing annotations...")

  # Link candidates to existing annotations via genomic overlap
  candidate_annotations <- link_candidates_to_annotations(con, candidate_import$file_id, background_file_id, verbose = verbose)

  if (verbose) {
    message("Import complete:")
    message("  - Candidate file ID: ", candidate_import$file_id)
    message("  - BED file created: ", length(candidate_bed), " genomic regions")
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
#' @export
extract_go_terms_for_enrichment <- function(con, foreground_file_id, background_file_id, 
                                           blast_param_id = NULL, verbose = TRUE) {

  if (verbose) message("Extracting GO terms for enrichment analysis...")

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
#' @export
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
  
  if (nrow(ontology_terms) == 0) {
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
  
  # Convert to funseqR format with additional clusterProfiler fields
  funseqr_results <- data.frame(
    go_id = cp_df$ID,
    go_term = cp_df$Description,
    go_category = ontology,
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
#' @export
extract_kegg_terms_for_enrichment <- function(con, foreground_file_id, background_file_id, 
                                            blast_param_id = NULL, verbose = TRUE) {

  if (verbose) message("Extracting KEGG pathways for enrichment analysis...")

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
#' @export
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
  
  if (nrow(pathways) == 0) {
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

#' Retrieve stored ORA results
#'
#' @param con Database connection object
#' @param analysis_id Integer. Analysis ID to retrieve
#' @param significance_filter Character. Filter by significance level. Default is NULL (no filter)
#'
#' @return List containing analysis metadata and results
#'
#' @export
get_ora_results <- function(con, analysis_id, significance_filter = NULL) {

  # Get analysis metadata
  analysis_query <- "SELECT * FROM ora_analyses WHERE analysis_id = ?"
  analysis_info <- DBI::dbGetQuery(con, analysis_query, list(analysis_id))

  if (nrow(analysis_info) == 0) {
    stop("No ORA analysis found with ID: ", analysis_id)
  }

  # Get results
  results_query <- "
    SELECT * FROM ora_results
    WHERE analysis_id = ?
    ORDER BY p_adjusted, fold_enrichment DESC
  "

  results <- DBI::dbGetQuery(con, results_query, list(analysis_id))

  # Apply significance filter if requested
  if (!is.null(significance_filter)) {
    results <- results[results$significance_level %in% significance_filter, ]
  }

  return(list(
    analysis_info = as.list(analysis_info[1, ]),
    results = results
  ))
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

    # Prepare results for insertion
    results_to_insert <- enrichment_results
    results_to_insert$analysis_id <- analysis_id

    # Determine if this is clusterProfiler results (has additional columns)
    has_clusterprofiler_cols <- all(c("gene_ratio", "bg_ratio", "qvalue", "gene_ids") %in% colnames(results_to_insert))
    
    # Map column names based on annotation type
    if (annotation_type == "GO") {
      term_id_col <- "go_id"
      term_name_col <- "go_term"
    } else {
      term_id_col <- "pathway_id" 
      term_name_col <- "pathway_name"
    }
    
    if (has_clusterprofiler_cols) {
      # clusterProfiler results with additional fields
      results_to_insert <- results_to_insert[, c("analysis_id", term_id_col, term_name_col,
                                                 "foreground_count", "background_count", "total_foreground",
                                                 "total_background", "expected_count", "fold_enrichment",
                                                 "p_value", "p_adjusted", "significance_level",
                                                 "gene_ratio", "bg_ratio", "qvalue", "gene_ids")]
      
      # Standardize column names for database
      colnames(results_to_insert)[2:3] <- c("term_id", "term_name")
      results_to_insert$annotation_type <- annotation_type
      results_to_insert$term_type <- term_type
      
      result_query <- "
        INSERT INTO ora_results
        (analysis_id, term_id, term_name, annotation_type, term_type, foreground_count, background_count,
         total_foreground, total_background, expected_count, fold_enrichment,
         p_value, p_adjusted, significance_level, gene_ratio, bg_ratio, qvalue, gene_ids)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
      "
    } else {
      # Legacy results
      results_to_insert <- results_to_insert[, c("analysis_id", term_id_col, term_name_col,
                                                 "foreground_count", "background_count", "total_foreground",
                                                 "total_background", "expected_count", "fold_enrichment",
                                                 "p_value", "p_adjusted", "significance_level")]
      
      # Standardize column names for database
      colnames(results_to_insert)[2:3] <- c("term_id", "term_name")
      results_to_insert$annotation_type <- annotation_type
      results_to_insert$term_type <- term_type
      
      result_query <- "
        INSERT INTO ora_results
        (analysis_id, term_id, term_name, annotation_type, term_type, foreground_count, background_count,
         total_foreground, total_background, expected_count, fold_enrichment,
         p_value, p_adjusted, significance_level)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
      "
    }

    for (i in 1:nrow(results_to_insert)) {
      row_data <- as.list(results_to_insert[i, ])
      names(row_data) <- NULL  # Remove names to use with anonymous placeholders
      DBI::dbExecute(con, result_query, row_data)
    }
  }

  if (verbose) {
    message("Storage complete:")
    message("  - Analysis ID: ", analysis_id)
    message("  - Results stored: ", nrow(enrichment_results))
  }

  return(analysis_id)
}

