#' Loci Traceability Functions for funseqR
#'
#' Functions to trace enriched GO terms and pathways back to specific genomic loci,
#' enabling downstream analysis such as environmental correlation mapping.

# INTERNAL HELPER FUNCTIONS (NOT EXPORTED)

#' Extract GO term IDs from enrichment results
#' @param enrichment_results Data frame from GO enrichment analysis
#' @param significance_threshold Numeric. FDR threshold for significance
#' @return Character vector of significant GO term IDs
#' @keywords internal
.extract_significant_go_ids <- function(enrichment_results, significance_threshold = 0.05) {
  if (is.null(enrichment_results) || nrow(enrichment_results) == 0) {
    return(character(0))
  }
  
  significant_terms <- enrichment_results[
    !is.na(enrichment_results$p_adjusted) & 
    enrichment_results$p_adjusted < significance_threshold, 
  ]
  
  return(unique(significant_terms$go_id))
}

#' Query database to trace GO terms to genomic loci
#' @param con Database connection
#' @param go_ids Character vector of GO term IDs
#' @param candidate_file_id Integer. File ID of candidate dataset
#' @return Data frame with genomic coordinates and functional annotations
#' @keywords internal
.query_go_terms_to_loci <- function(con, go_ids, candidate_file_id) {
  if (length(go_ids) == 0) {
    return(data.frame(
      chromosome = character(0),
      position = integer(0),
      ref = character(0),
      alt = character(0),
      go_id = character(0),
      go_term = character(0),
      go_category = character(0),
      gene_names = character(0),
      uniprot_accession = character(0),
      blast_evalue = numeric(0),
      blast_identity = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  
  # Create placeholders for SQL IN clause
  go_placeholders <- paste(rep("?", length(go_ids)), collapse = ",")
  
  query <- paste0("
    SELECT DISTINCT
      v.chromosome,
      v.position,
      v.ref,
      v.alt,
      gt.go_id,
      gt.go_term,
      gt.go_category,
      a.gene_names,
      a.uniprot_accession,
      br.e_value as blast_evalue,
      br.percent_identity as blast_identity
    FROM vcf_data v
    JOIN input_files if_cand ON v.file_id = if_cand.file_id
    JOIN flanking_sequences fs ON v.vcf_id = fs.vcf_id
    JOIN blast_results br ON fs.flanking_id = br.flanking_id
    JOIN annotations a ON br.blast_result_id = a.blast_result_id
    JOIN go_terms gt ON a.annotation_id = gt.annotation_id
    WHERE if_cand.file_id = ?
      AND gt.go_id IN (", go_placeholders, ")
    ORDER BY v.chromosome, v.position, gt.go_id
  ")
  
  params <- c(candidate_file_id, go_ids)
  result <- DBI::dbGetQuery(con, query, params = params)
  
  return(result)
}

#' Merge enrichment results with genomic coordinates
#' @param loci_data Data frame from database query
#' @param enrichment_results Data frame with enrichment statistics
#' @return Data frame with combined information
#' @keywords internal
.merge_loci_with_enrichment <- function(loci_data, enrichment_results) {
  if (nrow(loci_data) == 0 || nrow(enrichment_results) == 0) {
    return(loci_data)
  }
  
  # Merge by GO term ID to add enrichment statistics
  merged <- merge(
    loci_data,
    enrichment_results[, c("go_id", "fold_enrichment", "p_adjusted", "significance_level")],
    by = "go_id",
    all.x = TRUE
  )
  
  return(merged)
}

#' Create summary of loci per GO term
#' @param traced_data Data frame from trace analysis
#' @return Data frame summarizing loci counts per GO term
#' @keywords internal
.summarize_loci_per_term <- function(traced_data) {
  if (nrow(traced_data) == 0) {
    return(data.frame(
      go_id = character(0),
      go_term = character(0),
      go_category = character(0),
      loci_count = integer(0),
      chromosomes = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  summary_data <- aggregate(
    cbind(position, chromosome) ~ go_id + go_term + go_category,
    data = traced_data,
    FUN = function(x) length(unique(x))
  )
  
  names(summary_data)[names(summary_data) == "position"] <- "loci_count"
  
  # Add chromosome information
  chrom_summary <- aggregate(
    chromosome ~ go_id,
    data = traced_data,
    FUN = function(x) paste(sort(unique(x)), collapse = ", ")
  )
  
  summary_data <- merge(summary_data, chrom_summary, by = "go_id")
  summary_data$chromosome.y <- NULL
  names(summary_data)[names(summary_data) == "chromosome.y"] <- "chromosomes"
  
  # Reorder columns
  summary_data <- summary_data[, c("go_id", "go_term", "go_category", "loci_count", "chromosomes")]
  
  return(summary_data[order(summary_data$loci_count, decreasing = TRUE), ])
}

# EXPORTED FUNCTIONS

#' Trace enriched GO terms back to genomic loci
#'
#' This function traces significantly enriched GO terms from enrichment analysis
#' back to the specific genomic loci (chromosome positions) that contributed to
#' the enrichment signal. This enables downstream analysis such as environmental
#' correlation mapping and validation of functional predictions.
#'
#' @param con Database connection object
#' @param enrichment_results Data frame or list of enrichment results from GO analysis
#' @param candidate_file_id Integer. File ID of the candidate dataset used in enrichment analysis
#' @param significance_threshold Numeric. FDR threshold for selecting significant terms. Default is 0.05
#' @param include_blast_metrics Logical. Include BLAST quality metrics in output. Default is TRUE
#'
#' @return Data frame with columns:
#'   \itemize{
#'     \item chromosome: Chromosome name
#'     \item position: Genomic position
#'     \item ref: Reference allele
#'     \item alt: Alternative allele
#'     \item go_id: GO term identifier
#'     \item go_term: GO term description
#'     \item go_category: GO category (BP, MF, CC)
#'     \item gene_names: Associated gene names
#'     \item uniprot_accession: UniProt accession
#'     \item fold_enrichment: Enrichment fold change
#'     \item p_adjusted: FDR-adjusted p-value
#'     \item blast_evalue: BLAST e-value (if include_blast_metrics = TRUE)
#'     \item blast_identity: BLAST percent identity (if include_blast_metrics = TRUE)
#'   }
#'
#' @details
#' This function traverses the complete annotation pipeline in reverse:
#' GO terms → annotations → BLAST results → flanking sequences → VCF entries
#' 
#' It can handle enrichment results from individual ontologies or combined results
#' from multiple ontologies. The function automatically identifies significantly
#' enriched terms and traces them back to contributing genomic loci.
#'
#' @examples
#' \dontrun{
#' # Trace enrichment results back to genomic loci
#' con <- connect_funseq_db("analysis.db")
#' enrich_res <- run_go_enrichment_workflow(con, "candidates.vcf")
#' 
#' # Trace BP enrichment to loci
#' bp_loci <- trace_enriched_terms_to_loci(con, enrich_res$enrichment_results$BP, 
#'                                         enrich_res$candidate_import$file_id)
#' 
#' # Trace combined results
#' all_results <- do.call(rbind, enrich_res$enrichment_results)
#' all_loci <- trace_enriched_terms_to_loci(con, all_results, 
#'                                          enrich_res$candidate_import$file_id)
#' 
#' # Use for environmental correlation
#' write.table(all_loci[, c("chromosome", "position", "ref", "alt")], 
#'             "enriched_loci.txt", row.names = FALSE)
#' }
#'
#' @export
trace_enriched_terms_to_loci <- function(con, enrichment_results, candidate_file_id,
                                         significance_threshold = 0.05, 
                                         include_blast_metrics = TRUE) {
  
  # Validate inputs
  if (!is_db_connected(con)) {
    stop("Database connection is not valid")
  }
  
  if (is.null(enrichment_results) || length(enrichment_results) == 0) {
    stop("enrichment_results cannot be NULL or empty")
  }
  
  # Handle list input (e.g., from workflow results)
  if (is.list(enrichment_results) && !is.data.frame(enrichment_results)) {
    # Combine all enrichment results
    enrichment_results <- do.call(rbind, enrichment_results)
  }
  
  if (!is.data.frame(enrichment_results)) {
    stop("enrichment_results must be a data frame or list of data frames")
  }
  
  # Check required columns
  required_cols <- c("go_id", "go_term", "go_category", "p_adjusted")
  missing_cols <- setdiff(required_cols, names(enrichment_results))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in enrichment_results: ", paste(missing_cols, collapse = ", "))
  }
  
  # Extract significant GO term IDs
  significant_go_ids <- .extract_significant_go_ids(enrichment_results, significance_threshold)
  
  if (length(significant_go_ids) == 0) {
    warning("No significantly enriched terms found at threshold ", significance_threshold)
    return(data.frame())
  }
  
  # Query database to trace GO terms to loci
  loci_data <- .query_go_terms_to_loci(con, significant_go_ids, candidate_file_id)
  
  if (nrow(loci_data) == 0) {
    warning("No genomic loci found for enriched GO terms")
    return(data.frame())
  }
  
  # Merge with enrichment statistics
  result <- .merge_loci_with_enrichment(loci_data, enrichment_results)
  
  # Remove BLAST metrics if not requested
  if (!include_blast_metrics) {
    result$blast_evalue <- NULL
    result$blast_identity <- NULL
  }
  
  # Sort by significance and position
  result <- result[order(result$p_adjusted, result$chromosome, result$position), ]
  rownames(result) <- NULL
  
  return(result)
}

#' Export genomic coordinates of enriched loci
#'
#' Export genomic coordinates of loci associated with significantly enriched
#' functional terms in various formats suitable for downstream analysis.
#'
#' @param con Database connection object
#' @param enrichment_results Data frame or list of enrichment results
#' @param candidate_file_id Integer. File ID of the candidate dataset
#' @param output_file Character. Output file path (extension determines format)
#' @param format Character. Output format: "vcf", "bed", or "table". Default is "table"
#' @param significance_threshold Numeric. FDR threshold for significance. Default is 0.05
#' @param include_annotations Logical. Include functional annotations in output. Default is TRUE
#'
#' @return Invisible path to output file
#'
#' @details
#' Supported output formats:
#' \itemize{
#'   \item "table": Tab-delimited text file with all annotation information
#'   \item "bed": BED format file for use with genome browsers
#'   \item "vcf": VCF format file with enrichment information in INFO fields
#' }
#'
#' @examples
#' \dontrun{
#' # Export enriched loci coordinates
#' export_enriched_loci(con, enrichment_results, candidate_file_id, 
#'                      "enriched_loci.txt")
#' 
#' # Export as BED file for genome browser
#' export_enriched_loci(con, enrichment_results, candidate_file_id,
#'                      "enriched_loci.bed", format = "bed")
#' }
#'
#' @export
export_enriched_loci <- function(con, enrichment_results, candidate_file_id, output_file,
                                 format = "table", significance_threshold = 0.05,
                                 include_annotations = TRUE) {
  
  # Trace enriched terms to loci
  loci_data <- trace_enriched_terms_to_loci(
    con, enrichment_results, candidate_file_id,
    significance_threshold = significance_threshold,
    include_blast_metrics = include_annotations
  )
  
  if (nrow(loci_data) == 0) {
    stop("No enriched loci found to export")
  }
  
  # Auto-detect format from file extension if not specified
  if (missing(format)) {
    ext <- tools::file_ext(output_file)
    format <- switch(tolower(ext),
      "bed" = "bed",
      "vcf" = "vcf",
      "table"
    )
  }
  
  # Export based on format
  switch(format,
    "table" = {
      write.table(loci_data, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
    },
    "bed" = {
      bed_data <- data.frame(
        chrom = loci_data$chromosome,
        chromStart = loci_data$position - 1,  # BED is 0-based
        chromEnd = loci_data$position,
        name = paste0(loci_data$go_id, ":", gsub(" ", "_", loci_data$go_term)),
        score = round(-log10(loci_data$p_adjusted) * 100),
        strand = "."
      )
      write.table(bed_data, output_file, sep = "\t", row.names = FALSE, 
                  col.names = FALSE, quote = FALSE)
    },
    "vcf" = {
      # Create basic VCF structure
      vcf_data <- data.frame(
        CHROM = loci_data$chromosome,
        POS = loci_data$position,
        ID = ".",
        REF = loci_data$ref,
        ALT = loci_data$alt,
        QUAL = ".",
        FILTER = "PASS",
        INFO = paste0("GO_ID=", loci_data$go_id, 
                     ";GO_TERM=", gsub(" ", "_", loci_data$go_term),
                     ";GO_CAT=", loci_data$go_category,
                     ";FOLD_ENRICH=", round(loci_data$fold_enrichment, 2),
                     ";FDR=", format(loci_data$p_adjusted, scientific = TRUE))
      )
      
      # Write VCF header
      cat("##fileformat=VCFv4.2\n", file = output_file)
      cat("##source=funseqR\n", file = output_file, append = TRUE)
      cat("##INFO=<ID=GO_ID,Number=1,Type=String,Description=\"GO term ID\">\n", 
          file = output_file, append = TRUE)
      cat("##INFO=<ID=GO_TERM,Number=1,Type=String,Description=\"GO term description\">\n", 
          file = output_file, append = TRUE)
      cat("##INFO=<ID=GO_CAT,Number=1,Type=String,Description=\"GO category\">\n", 
          file = output_file, append = TRUE)
      cat("##INFO=<ID=FOLD_ENRICH,Number=1,Type=Float,Description=\"Fold enrichment\">\n", 
          file = output_file, append = TRUE)
      cat("##INFO=<ID=FDR,Number=1,Type=Float,Description=\"FDR-adjusted p-value\">\n", 
          file = output_file, append = TRUE)
      
      # Write data
      write.table(vcf_data, output_file, sep = "\t", row.names = FALSE, quote = FALSE, append = TRUE)
    },
    stop("Unsupported format: ", format, ". Use 'table', 'bed', or 'vcf'")
  )
  
  message("Exported ", nrow(loci_data), " enriched loci to: ", output_file)
  invisible(output_file)
}

#' Summarize functional annotations at genomic loci
#'
#' Create comprehensive summaries of functional annotations associated with
#' enriched genomic loci, including statistics about GO terms, pathways, and
#' genomic distribution.
#'
#' @param con Database connection object
#' @param enrichment_results Data frame or list of enrichment results
#' @param candidate_file_id Integer. File ID of the candidate dataset
#' @param significance_threshold Numeric. FDR threshold for significance. Default is 0.05
#'
#' @return List containing:
#'   \itemize{
#'     \item loci_summary: Data frame with per-locus functional annotations
#'     \item term_summary: Data frame with per-term genomic distribution
#'     \item chromosome_summary: Data frame with per-chromosome enrichment counts
#'     \item overview: Named vector with summary statistics
#'   }
#'
#' @examples
#' \dontrun{
#' # Create functional summary of enriched loci
#' summary_results <- summarize_functional_loci(con, enrichment_results, candidate_file_id)
#' 
#' # View overview statistics
#' print(summary_results$overview)
#' 
#' # View per-chromosome distribution
#' print(summary_results$chromosome_summary)
#' }
#'
#' @export
summarize_functional_loci <- function(con, enrichment_results, candidate_file_id,
                                     significance_threshold = 0.05) {
  
  # Trace enriched terms to loci
  loci_data <- trace_enriched_terms_to_loci(
    con, enrichment_results, candidate_file_id,
    significance_threshold = significance_threshold,
    include_blast_metrics = TRUE
  )
  
  if (nrow(loci_data) == 0) {
    return(list(
      loci_summary = data.frame(),
      term_summary = data.frame(),
      chromosome_summary = data.frame(),
      overview = c(total_loci = 0, total_terms = 0, total_chromosomes = 0)
    ))
  }
  
  # Per-locus summary (one row per unique genomic position)
  loci_summary <- aggregate(
    cbind(go_id, go_term, go_category) ~ chromosome + position + ref + alt + gene_names + uniprot_accession,
    data = loci_data,
    FUN = function(x) length(unique(x))
  )
  names(loci_summary)[names(loci_summary) %in% c("go_id", "go_term", "go_category")] <- 
    c("num_go_terms", "num_unique_terms", "num_categories")
  
  # Add functional term lists
  term_lists <- aggregate(
    go_term ~ chromosome + position + ref + alt,
    data = loci_data,
    FUN = function(x) paste(unique(x), collapse = "; ")
  )
  loci_summary <- merge(loci_summary, term_lists, by = c("chromosome", "position", "ref", "alt"))
  names(loci_summary)[names(loci_summary) == "go_term"] <- "enriched_terms"
  
  # Per-term summary (one row per GO term)
  term_summary <- .summarize_loci_per_term(loci_data)
  
  # Per-chromosome summary
  chromosome_summary <- aggregate(
    cbind(position, go_id) ~ chromosome,
    data = loci_data,
    FUN = function(x) length(unique(x))
  )
  names(chromosome_summary) <- c("chromosome", "unique_loci", "enriched_terms")
  
  # Add GO category breakdown per chromosome
  category_breakdown <- aggregate(
    go_category ~ chromosome,
    data = loci_data,
    FUN = function(x) {
      counts <- table(x)
      paste(paste(names(counts), counts, sep = ":"), collapse = ", ")
    }
  )
  chromosome_summary <- merge(chromosome_summary, category_breakdown, by = "chromosome")
  names(chromosome_summary)[names(chromosome_summary) == "go_category"] <- "category_breakdown"
  
  # Overview statistics
  overview <- c(
    total_loci = length(unique(paste(loci_data$chromosome, loci_data$position))),
    total_terms = length(unique(loci_data$go_id)),
    total_chromosomes = length(unique(loci_data$chromosome)),
    bp_terms = sum(loci_data$go_category == "BP"),
    mf_terms = sum(loci_data$go_category == "MF"),
    cc_terms = sum(loci_data$go_category == "CC")
  )
  
  return(list(
    loci_summary = loci_summary[order(loci_summary$chromosome, loci_summary$position), ],
    term_summary = term_summary,
    chromosome_summary = chromosome_summary[order(chromosome_summary$unique_loci, decreasing = TRUE), ],
    overview = overview
  ))
}