# GO Enrichment Analysis Example
# 
# This example demonstrates how to use the new GO enrichment functionality
# to compare candidate adaptive loci vs background SNPs.
#
# Prerequisites:
# 1. A completed funseqR analysis with BLAST and annotations
# 2. Two VCF files: full dataset (background) and candidate subset (foreground)

library(funseqR)

# Example workflow for the snapper project
# =======================================

# Path setup
project_dir <- "/Users/brau0037/Library/CloudStorage/GoogleDrive-pygmyperch@gmail.com/My Drive/projects/snapperSG/south_coast/annotation/funseq_annotationMS"
db_path <- file.path(project_dir, "funseq_project.db")
candidate_vcf <- file.path(project_dir, "SA448_855.vcf")  # 855 candidate SNPs
background_vcf_id <- 1  # File ID of SA448_14699.vcf (full dataset)

# Connect to database
con <- connect_funseq_db(db_path)

# Method 1: Complete Workflow (Recommended)
# ==========================================

# Run the complete workflow - this handles everything automatically
results <- run_go_enrichment_workflow(
  con = con,
  project_id = 1,
  candidate_vcf_file = candidate_vcf,
  background_file_id = background_vcf_id,
  ontologies = c("BP", "MF"),  # Test Biological Process and Molecular Function
  min_genes = 5,               # Minimum genes required for testing a GO term
  max_genes = 500,             # Maximum genes (excludes very broad terms)
  store_results = TRUE,        # Store in database for future retrieval
  create_plots = TRUE,         # Generate visualization plots
  verbose = TRUE
)

# View summary
print(results$summary)

# Display plots
if (!is.null(results$plots$BP_bubble)) {
  print(results$plots$BP_bubble)
}

if (!is.null(results$plots$MF_bubble)) {
  print(results$plots$MF_bubble)
}

# Show summary table
if (!is.null(results$plots$BP_table)) {
  print(results$plots$BP_table)
}

# Method 2: Step-by-Step Workflow
# ================================

# If you want more control, you can run each step manually:

# Step 1: Import candidate file and link to annotations
candidate_import <- import_candidate_loci(con, 1, candidate_vcf, background_vcf_id)

# Step 2: Extract GO terms for both datasets
go_data <- extract_go_terms_for_enrichment(con, candidate_import$file_id, background_vcf_id)

# Step 3: Perform enrichment analysis
bp_results <- perform_go_enrichment(go_data, "BP", min_genes = 5)
mf_results <- perform_go_enrichment(go_data, "MF", min_genes = 5)

# Step 4: Create visualizations
bp_bubble <- create_go_bubble_plot(bp_results, max_terms = 20)
mf_bubble <- create_go_bubble_plot(mf_results, max_terms = 15)

# Optional: Create treemap (requires treemapify package)
if (requireNamespace("treemapify", quietly = TRUE)) {
  bp_treemap <- create_go_treemap(bp_results)
  print(bp_treemap)
}

# Step 5: Store results in database
bp_enrichment_id <- store_go_enrichment_results(con, 1, candidate_import$file_id, 
                                               background_vcf_id, bp_results, "BP")

mf_enrichment_id <- store_go_enrichment_results(con, 1, candidate_import$file_id, 
                                               background_vcf_id, mf_results, "MF")

# Retrieving Stored Results
# =========================

# You can retrieve stored results later:
stored_bp <- get_go_enrichment_results(con, bp_enrichment_id)
stored_mf <- get_go_enrichment_results(con, mf_enrichment_id)

# View analysis metadata
print(stored_bp$analysis_info)

# View significant results only
significant_bp <- get_go_enrichment_results(con, bp_enrichment_id, 
                                           significance_filter = c("significant", "highly_significant"))

# Creating Publication-Ready Tables
# ==================================

# Create formatted summary table
summary_table <- create_go_summary_table(bp_results, max_terms = 10)
print(summary_table)

# Create comparison plot across ontologies
comparison_plot <- create_go_comparison_plot(bp_results, mf_results, top_n = 8)
print(comparison_plot)

# Export results to files
# =======================

# Save plots
ggsave("bp_enrichment_bubble.png", bp_bubble, width = 12, height = 8, dpi = 300)
ggsave("mf_enrichment_bubble.png", mf_bubble, width = 12, height = 6, dpi = 300)

# Save summary table as CSV
write.csv(summary_table, "go_enrichment_summary.csv", row.names = FALSE)

# Save detailed results
write.csv(bp_results, "bp_enrichment_detailed.csv", row.names = FALSE)
write.csv(mf_results, "mf_enrichment_detailed.csv", row.names = FALSE)

# Interactive Analysis (if plotly available)
# ==========================================

if (requireNamespace("plotly", quietly = TRUE)) {
  interactive_plot <- create_interactive_go_plot(bp_results, max_terms = 25)
  # Save as HTML
  htmlwidgets::saveWidget(interactive_plot, "interactive_go_enrichment.html")
}

# Quality Control and Interpretation
# ===================================

cat("\n=== Analysis Summary ===\n")
cat("Candidate genes with GO annotations:", length(go_data$foreground$genes), "\n")
cat("Background genes with GO annotations:", length(go_data$background$genes), "\n")
cat("Biological Process terms tested:", nrow(bp_results), "\n")
cat("Biological Process terms significant (FDR < 0.05):", sum(bp_results$p_adjusted < 0.05), "\n")
cat("Molecular Function terms tested:", nrow(mf_results), "\n")
cat("Molecular Function terms significant (FDR < 0.05):", sum(mf_results$p_adjusted < 0.05), "\n")

# Show top 5 most enriched terms
cat("\n=== Top 5 Enriched Biological Processes ===\n")
top_bp <- bp_results[bp_results$p_adjusted < 0.05, ][1:min(5, sum(bp_results$p_adjusted < 0.05)), ]
for (i in 1:nrow(top_bp)) {
  cat(sprintf("%d. %s (Fold: %.2f, FDR: %.2e)\n", 
              i, top_bp$go_term[i], top_bp$fold_enrichment[i], top_bp$p_adjusted[i]))
}

cat("\n=== Top 5 Enriched Molecular Functions ===\n")
top_mf <- mf_results[mf_results$p_adjusted < 0.05, ][1:min(5, sum(mf_results$p_adjusted < 0.05)), ]
for (i in 1:nrow(top_mf)) {
  cat(sprintf("%d. %s (Fold: %.2f, FDR: %.2e)\n", 
              i, top_mf$go_term[i], top_mf$fold_enrichment[i], top_mf$p_adjusted[i]))
}

# Close connection
close_funseq_db(con)

cat("\nGO enrichment analysis complete!\n")
cat("Check the generated plots and tables for biological insights.\n")