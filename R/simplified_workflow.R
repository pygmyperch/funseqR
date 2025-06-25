#' Simplified funseqR Workflow Documentation
#'
#' This file documents the new simplified workflow for functional annotation analysis.
#'
#' @section New Simplified Workflow:
#'
#' The funseqR package now provides a streamlined workflow that reduces complexity
#' while maintaining all analytical capabilities:
#'
#' \subsection{Step 1: Run BLAST and Annotation (unchanged)}{
#' ```r
#' # Run BLAST analysis (same as before)
#' blast_results <- run_blast(sequences, database)
#' 
#' # Annotate BLAST results (same as before)
#' con <- annotate_blast_results(blast_results)
#' ```
#' }
#'
#' \subsection{Step 2: Process Annotations (NEW - replaces multiple functions)}{
#' ```r
#' # Process all functional annotations into standardized table
#' annotations <- process_annotations(
#'   con, 
#'   include = c("GO", "KEGG"),           # Default: both
#'   export_csv = "my_annotations.csv"   # Optional CSV export
#' )
#' 
#' # The result is a comprehensive data frame with:
#' # - locus_id, chromosome, position
#' # - gene_name, protein_name, uniprot_accession
#' # - GO terms, names, and categories
#' # - KEGG pathways with enhanced names
#' # - KEGG BRITE hierarchy classifications
#' # - Functional module assignments
#' ```
#' }
#'
#' \subsection{Step 3: Run Enrichment Analysis (simplified)}{
#' ```r
#' # Define candidate loci
#' candidates <- data.frame(
#'   chromosome = c("LG1", "LG2"), 
#'   position = c(12345, 67890)
#' )
#' 
#' # GO enrichment analysis
#' go_results <- run_go_enrichment_analysis(
#'   annotations, 
#'   candidate_loci = candidates,
#'   ontologies = c("BP", "MF", "CC")
#' )
#' 
#' # KEGG pathway enrichment analysis
#' kegg_results <- run_kegg_enrichment_analysis(
#'   annotations,
#'   candidate_loci = candidates
#' )
#' ```
#' }
#'
#' @section Key Improvements:
#' 
#' \itemize{
#'   \item \strong{Single function replaces many}: `process_annotations()` replaces 
#'     `summarize_go_annotations()`, `summarize_kegg_pathways()`, 
#'     `generate_candidate_loci_profile()`, and others
#'   \item \strong{Standardized output}: One consistent data frame format for all downstream analyses
#'   \item \strong{Always enhanced}: KEGGREST and ggkegg enhancements are built-in, not optional
#'   \item \strong{User-friendly}: CSV export for external analysis tools
#'   \item \strong{Future-proof}: Easy to add new annotation types (eggNOG, COG, etc.)
#'   \item \strong{Maintains power}: All statistical capabilities preserved
#' }
#'
#' @section Migration from Old Workflow:
#'
#' \subsection{Before (complex):}{
#' ```r
#' # Old workflow required many functions
#' go_summary <- summarize_go_annotations(con, candidate_loci = candidates)
#' kegg_summary <- summarize_kegg_pathways(con, candidate_loci = candidates, 
#'                                        use_keggrest_enhancement = TRUE)
#' profile <- generate_candidate_loci_profile(con, candidate_loci = candidates)
#' go_enrichment <- extract_go_terms_for_enrichment(con, ...)
#' # ... many more steps
#' ```
#' }
#'
#' \subsection{After (simple):}{
#' ```r
#' # New workflow: just 3 main steps
#' annotations <- process_annotations(con, include = c("GO", "KEGG"))
#' go_results <- run_go_enrichment_analysis(annotations, candidates)
#' kegg_results <- run_kegg_enrichment_analysis(annotations, candidates)
#' ```
#' }
#'
#' @section Output Format:
#'
#' The `process_annotations()` function returns a standardized data frame with these columns:
#' 
#' \describe{
#'   \item{locus_id}{Unique identifier: "file_chromosome_position"}
#'   \item{chromosome}{Chromosome name}
#'   \item{position}{Genomic position}
#'   \item{gene_name}{Primary gene name from UniProt}
#'   \item{protein_name}{Protein description}
#'   \item{uniprot_accession}{UniProt accession number(s)}
#'   \item{go_terms}{Semicolon-separated GO term IDs}
#'   \item{go_names}{Semicolon-separated GO term names}
#'   \item{go_categories}{Semicolon-separated GO categories (BP/MF/CC)}
#'   \item{kegg_pathways}{Semicolon-separated KEGG pathway IDs}
#'   \item{kegg_pathway_names}{Semicolon-separated enhanced pathway names}
#'   \item{kegg_brite_categories}{Semicolon-separated KEGG BRITE classifications}
#'   \item{kegg_modules}{Semicolon-separated functional module assignments}
#' }
#'
#' This format enables easy filtering, analysis, and export to other tools.
#'
#' @name simplified_workflow
NULL