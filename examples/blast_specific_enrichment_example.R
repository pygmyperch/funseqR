#' BLAST-Specific GO Enrichment Analysis Example
#' 
#' This example demonstrates how to compare functional enrichment results
#' from different BLAST runs (e.g., ORF sequences vs raw sequences) using
#' the simplified funseqR workflow with single blast_param_id parameter.

library(funseqR)

# Connect to database
con <- connect_funseq_db("analysis.db")

# Assume you have run multiple BLAST analyses
# Example scenario from user's use case:
# blast_results1 <- perform_blast_db(con, vcf_import$file_id, db_path, db_name,
#                                     blast_type = "diamond_blastx", seq_type = "orf_nuc", ...)
# blast_results2 <- perform_blast_db(con, vcf_import$file_id, db_path, db_name,
#                                     blast_type = "diamond_blastx", seq_type = "raw", ...)

# Check available BLAST runs
cat("=== Available BLAST Runs ===\n")
blast_runs <- DBI::dbGetQuery(con, "
  SELECT blast_param_id, blast_type, db_name, execution_date 
  FROM blast_parameters 
  ORDER BY execution_date DESC
")
print(blast_runs)

# Count annotations per BLAST run  
cat("\n=== Annotations per BLAST Run ===\n")
annotation_counts <- DBI::dbGetQuery(con, "
  SELECT bp.blast_param_id, bp.db_name, bp.execution_date, COUNT(a.annotation_id) as annotation_count
  FROM blast_parameters bp
  LEFT JOIN blast_results br ON bp.blast_param_id = br.blast_param_id  
  LEFT JOIN annotations a ON br.blast_result_id = a.blast_result_id
  GROUP BY bp.blast_param_id
  ORDER BY bp.execution_date DESC
")
print(annotation_counts)

# === COMPARING ORF vs RAW SEQUENCE ENRICHMENT ===

# Assuming BLAST run IDs from the user's scenario:
orf_blast_id <- 1    # ORF sequences vs teleost proteins
raw_blast_id <- 2    # Raw sequences vs teleost proteins
candidate_vcf_file <- "candidates.vcf"

cat("\n=== Running ORF-based GO Enrichment ===\n")
# Run enrichment using only ORF-based annotations
enrich_res_orf <- run_go_enrichment_workflow(
  con, 
  candidate_vcf_file, 
  blast_param_id = orf_blast_id,  # Use ORF-based BLAST run for both candidate and background
  ontologies = c("BP", "MF", "CC"),
  verbose = TRUE
)

cat("\n=== Running Raw Sequence GO Enrichment ===\n")
# Run enrichment using only raw sequence annotations
enrich_res_raw <- run_go_enrichment_workflow(
  con, 
  candidate_vcf_file,
  blast_param_id = raw_blast_id,  # Use raw sequence BLAST run for both candidate and background
  ontologies = c("BP", "MF", "CC"),
  verbose = TRUE
)

cat("\n=== Running Combined Enrichment (All Annotations) ===\n")
# Run enrichment using all available annotations (original behavior)
enrich_res_all <- run_go_enrichment_workflow(
  con, 
  candidate_vcf_file,
  # No blast_param_id specified = use all annotations
  ontologies = c("BP", "MF", "CC"),
  verbose = TRUE
)

# === COMPARE RESULTS ===

cat("\n=== Comparison Summary ===\n")
cat("ORF-based enrichment:\n")
cat("  - Significant BP terms:", sum(enrich_res_orf$enrichment_results$BP$p_adjusted < 0.05, na.rm = TRUE), "\n")
cat("  - Significant MF terms:", sum(enrich_res_orf$enrichment_results$MF$p_adjusted < 0.05, na.rm = TRUE), "\n")
cat("  - Significant CC terms:", sum(enrich_res_orf$enrichment_results$CC$p_adjusted < 0.05, na.rm = TRUE), "\n")

cat("\nRaw sequence enrichment:\n")
cat("  - Significant BP terms:", sum(enrich_res_raw$enrichment_results$BP$p_adjusted < 0.05, na.rm = TRUE), "\n")
cat("  - Significant MF terms:", sum(enrich_res_raw$enrichment_results$MF$p_adjusted < 0.05, na.rm = TRUE), "\n")
cat("  - Significant CC terms:", sum(enrich_res_raw$enrichment_results$CC$p_adjusted < 0.05, na.rm = TRUE), "\n")

cat("\nCombined enrichment:\n")
cat("  - Significant BP terms:", sum(enrich_res_all$enrichment_results$BP$p_adjusted < 0.05, na.rm = TRUE), "\n")
cat("  - Significant MF terms:", sum(enrich_res_all$enrichment_results$MF$p_adjusted < 0.05, na.rm = TRUE), "\n")
cat("  - Significant CC terms:", sum(enrich_res_all$enrichment_results$CC$p_adjusted < 0.05, na.rm = TRUE), "\n")

# === CREATE SUMMARY TABLES ===

# ORF-based enrichment tables
mf_table_orf <- create_go_summary_table(enrich_res_orf$enrichment_results$MF)
cc_table_orf <- create_go_summary_table(enrich_res_orf$enrichment_results$CC)

# Raw sequence enrichment tables  
mf_table_raw <- create_go_summary_table(enrich_res_raw$enrichment_results$MF)
cc_table_raw <- create_go_summary_table(enrich_res_raw$enrichment_results$CC)

cat("\n=== Top ORF-based MF Terms ===\n")
print(head(mf_table_orf))

cat("\n=== Top Raw Sequence MF Terms ===\n") 
print(head(mf_table_raw))

# === TRACEABILITY ANALYSIS WITH BLAST FILTERING ===

cat("\n=== Tracing ORF-based Enrichment to Genomic Loci ===\n")
# Trace ORF enrichment back to specific genomic loci
orf_loci <- trace_enriched_terms_to_loci(
  con, 
  enrichment_results = enrich_res_orf$enrichment_results,
  candidate_file_id = enrich_res_orf$candidate_import$file_id,
  blast_param_id = orf_blast_id  # Only trace ORF-based annotations
)

cat("ORF enrichment traces to", nrow(orf_loci), "genomic loci\n")

cat("\n=== Tracing Raw Sequence Enrichment to Genomic Loci ===\n")
# Trace raw sequence enrichment back to specific genomic loci  
raw_loci <- trace_enriched_terms_to_loci(
  con,
  enrichment_results = enrich_res_raw$enrichment_results,
  candidate_file_id = enrich_res_raw$candidate_import$file_id,
  blast_param_id = raw_blast_id  # Only trace raw sequence annotations
)

cat("Raw sequence enrichment traces to", nrow(raw_loci), "genomic loci\n")

# === FUNCTIONAL SUMMARIES ===

cat("\n=== ORF-based Functional Summary ===\n")
orf_summary <- summarize_functional_loci(
  con,
  enrichment_results = enrich_res_orf$enrichment_results,
  candidate_file_id = enrich_res_orf$candidate_import$file_id,
  blast_param_id = orf_blast_id
)
print(orf_summary$overview)

cat("\n=== Raw Sequence Functional Summary ===\n")
raw_summary <- summarize_functional_loci(
  con,
  enrichment_results = enrich_res_raw$enrichment_results, 
  candidate_file_id = enrich_res_raw$candidate_import$file_id,
  blast_param_id = raw_blast_id
)
print(raw_summary$overview)

# === EXPORT BLAST-SPECIFIC LOCI ===

# Export ORF-based enriched loci for environmental analysis
export_enriched_loci(
  con,
  enrichment_results = enrich_res_orf$enrichment_results,
  candidate_file_id = enrich_res_orf$candidate_import$file_id,
  output_file = "orf_enriched_loci.txt",
  blast_param_id = orf_blast_id,
  format = "table"
)

# Export raw sequence enriched loci
export_enriched_loci(
  con,
  enrichment_results = enrich_res_raw$enrichment_results,
  candidate_file_id = enrich_res_raw$candidate_import$file_id,
  output_file = "raw_enriched_loci.txt", 
  blast_param_id = raw_blast_id,
  format = "table"
)

# === VISUALIZATION COMPARISON ===

# Create comparison plots
if (requireNamespace("ggplot2", quietly = TRUE)) {
  # ORF-based bubble plot
  orf_plot <- create_go_bubble_plot(enrich_res_orf$enrichment_results$MF, 
                                    title = "ORF-based MF Enrichment")
  print(orf_plot)
  ggsave("orf_mf_enrichment.png", orf_plot, width = 12, height = 8)
  
  # Raw sequence bubble plot
  raw_plot <- create_go_bubble_plot(enrich_res_raw$enrichment_results$MF,
                                    title = "Raw Sequence MF Enrichment") 
  print(raw_plot)
  ggsave("raw_mf_enrichment.png", raw_plot, width = 12, height = 8)
}

# === MANHATTAN PLOTS WITH BLAST FILTERING ===

# Example GEA results (replace with your actual environmental analysis results)
gea_results <- data.frame(
  chromosome = c("1", "1", "2", "3", "1", "2"),
  position = c(1234567, 2345678, 1567890, 890123, 3456789, 2123456),
  p_value = c(0.001, 0.0001, 0.005, 0.01, 0.0005, 0.002),
  effect_size = c(0.8, 1.2, 0.6, 0.9, 1.1, 0.7)
)

# Create Manhattan plots for each annotation strategy
if (requireNamespace("ggplot2", quietly = TRUE)) {
  # ORF-based functional Manhattan plot
  # Note: You would need to trace the enrichment results and merge with GEA data
  # This is a conceptual example showing the workflow
  
  cat("\n=== Creating ORF-based Manhattan Plot ===\n")
  # manhattan_orf <- create_functional_manhattan_plot(
  #   con, enrich_res_orf$enrichment_results,
  #   enrich_res_orf$candidate_import$file_id,
  #   gea_results, y_metric = "p_value",
  #   title = "ORF-based Functional Manhattan Plot"
  # )
  
  cat("\n=== Creating Raw Sequence Manhattan Plot ===\n") 
  # manhattan_raw <- create_functional_manhattan_plot(
  #   con, enrich_res_raw$enrichment_results,
  #   enrich_res_raw$candidate_import$file_id,
  #   gea_results, y_metric = "p_value", 
  #   title = "Raw Sequence Functional Manhattan Plot"
  # )
}

cat("\n=== Analysis Complete ===\n")
cat("Files generated:\n")
cat("  - orf_enriched_loci.txt: Genomic coordinates from ORF-based enrichment\n")
cat("  - raw_enriched_loci.txt: Genomic coordinates from raw sequence enrichment\n")
cat("  - orf_mf_enrichment.png: ORF-based MF enrichment bubble plot\n")
cat("  - raw_mf_enrichment.png: Raw sequence MF enrichment bubble plot\n")

# Close database connection
close_funseq_db(con)

cat("\nThis analysis demonstrates how to:\n")
cat("1. Compare functional profiles from different sequence types\n")
cat("2. Trace enrichments back to specific BLAST runs\n") 
cat("3. Export annotation-specific loci for environmental analysis\n")
cat("4. Validate that enriched functions come from expected annotation sources\n")