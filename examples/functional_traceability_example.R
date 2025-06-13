#' Functional Traceability Analysis Example
#' 
#' This example demonstrates how to trace enriched GO terms back to genomic loci
#' and create publication-ready visualizations linking functional predictions
#' to environmental genomic analysis results.

library(funseqR)

# Connect to database
con <- connect_funseq_db("analysis.db")

# Assume you have already run GO enrichment workflow
# enrich_res <- run_go_enrichment_workflow(con, candidate_vcf_file, 
#                                          ontologies = c("BP", "MF", "CC"))

# For this example, we'll assume you have enrichment results like:
# mf_table <- create_go_summary_table(enrich_res$enrichment_results$MF)
# cc_table <- create_go_summary_table(enrich_res$enrichment_results$CC)

# candidate_file_id <- enrich_res$candidate_import$file_id

# Example enrichment results and file ID (replace with your actual data)
candidate_file_id <- 1  # Your actual candidate file ID

# === STEP 1: Trace enriched terms back to genomic loci ===

# Trace all significant enriched terms to their genomic origins
all_enriched_loci <- trace_enriched_terms_to_loci(
  con, 
  enrichment_results = enrich_res$enrichment_results,  # Use your actual results
  candidate_file_id = candidate_file_id,
  significance_threshold = 0.05,
  include_blast_metrics = TRUE
)

print("Genomic loci driving functional enrichment:")
print(head(all_enriched_loci))

# === STEP 2: Export enriched loci for environmental analysis ===

# Export as simple table for correlation analysis
export_enriched_loci(
  con, 
  enrichment_results = enrich_res$enrichment_results,
  candidate_file_id = candidate_file_id,
  output_file = "enriched_loci_for_environmental_analysis.txt",
  format = "table"
)

# Export as BED file for genome browser visualization
export_enriched_loci(
  con, 
  enrichment_results = enrich_res$enrichment_results,
  candidate_file_id = candidate_file_id,
  output_file = "enriched_loci.bed",
  format = "bed"
)

# === STEP 3: Create functional summary ===

functional_summary <- summarize_functional_loci(
  con,
  enrichment_results = enrich_res$enrichment_results,
  candidate_file_id = candidate_file_id
)

print("Overview of functional enrichment:")
print(functional_summary$overview)

print("Per-chromosome distribution:")
print(functional_summary$chromosome_summary)

print("Most functionally enriched loci:")
print(head(functional_summary$loci_summary))

# === STEP 4: Create functional Manhattan plot ===

# Example GEA results (replace with your actual environmental analysis results)
gea_results <- data.frame(
  chromosome = c("1", "1", "2", "3", "1", "2"),
  position = c(1234567, 2345678, 1567890, 890123, 3456789, 2123456),
  p_value = c(0.001, 0.0001, 0.005, 0.01, 0.0005, 0.002),
  effect_size = c(0.8, 1.2, 0.6, 0.9, 1.1, 0.7)
)

# Create Manhattan plot linking functional enrichment to environmental signals
manhattan_plot <- create_functional_manhattan_plot(
  con,
  enrichment_results = enrich_res$enrichment_results,
  candidate_file_id = candidate_file_id,
  gea_results = gea_results,
  y_metric = "p_value",
  y_transform = "-log10",
  color_by = "go_category",
  point_size_by = "fold_enrichment",
  max_labels = 15,
  title = "Environmental Association vs Functional Enrichment"
)

print(manhattan_plot)

# Save the plot
ggsave("functional_manhattan_plot.png", manhattan_plot, 
       width = 14, height = 8, dpi = 300)

# === STEP 5: Validate functional predictions ===

# Check if functionally enriched loci are also environmentally significant
validation_analysis <- merge(
  all_enriched_loci[, c("chromosome", "position", "go_term", "fold_enrichment", "p_adjusted")],
  gea_results,
  by = c("chromosome", "position"),
  all.x = TRUE
)

# Correlation between functional enrichment and environmental significance
correlation_result <- cor.test(
  -log10(validation_analysis$p_adjusted), 
  -log10(validation_analysis$p_value),
  use = "complete.obs"
)

print(paste("Correlation between functional and environmental significance:", 
            round(correlation_result$estimate, 3)))
print(paste("P-value for correlation:", format(correlation_result$p.value, scientific = TRUE)))

# === STEP 6: Focus on specific functional categories ===

# Extract loci associated with ion channel functions (example)
channel_terms <- all_enriched_loci[
  grepl("channel|transport", all_enriched_loci$go_term, ignore.case = TRUE),
]

if (nrow(channel_terms) > 0) {
  print("Ion channel and transport-related loci:")
  print(channel_terms[, c("chromosome", "position", "go_term", "fold_enrichment")])
  
  # Export these specific loci for targeted environmental analysis
  write.table(
    channel_terms[, c("chromosome", "position", "ref", "alt", "go_term")],
    "ion_channel_loci.txt",
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}

# === STEP 7: Regional analysis of enriched regions ===

# If you want to examine a specific genomic region in detail
if (nrow(all_enriched_loci) > 0) {
  # Example: zoom into first enriched chromosome
  target_chrom <- all_enriched_loci$chromosome[1]
  target_positions <- all_enriched_loci[all_enriched_loci$chromosome == target_chrom, "position"]
  
  # Define region around enriched loci (±500kb)
  region_start <- min(target_positions) - 500000
  region_end <- max(target_positions) + 500000
  
  regional_plot <- create_enrichment_region_plot(
    con,
    enrichment_results = enrich_res$enrichment_results,
    candidate_file_id = candidate_file_id,
    chromosome = target_chrom,
    start = max(1, region_start),  # Ensure positive coordinates
    end = region_end
  )
  
  print(regional_plot)
  ggsave(paste0("enrichment_region_", target_chrom, ".png"), regional_plot,
         width = 12, height = 6, dpi = 300)
}

# === STEP 8: Summary for publication ===

cat("\n=== FUNCTIONAL TRACEABILITY SUMMARY ===\n")
cat("Total enriched loci identified:", nrow(all_enriched_loci), "\n")
cat("Chromosomes with enriched loci:", length(unique(all_enriched_loci$chromosome)), "\n")
cat("Unique GO terms represented:", length(unique(all_enriched_loci$go_id)), "\n")
cat("Files exported for environmental analysis:\n")
cat("  - enriched_loci_for_environmental_analysis.txt\n")
cat("  - enriched_loci.bed\n")
cat("  - functional_manhattan_plot.png\n")

# Close database connection
close_funseq_db(con)