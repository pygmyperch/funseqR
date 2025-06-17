#' Manhattan Plot Example with Functional Annotation
#'
#' This example demonstrates the updated Manhattan plot functionality
#' with chromosome consolidation and functional enrichment highlighting.

library(funseqR)

# Connect to database
con <- connect_funseq_db("funseq_project.db")

# Load your RDA results (replace with your actual data loading)
# load("path/to/your/rda_results.RData")  # Should contain rda.simple.pq

# === EXAMPLE USAGE ===

# 1. Create basic Manhattan plot without functional annotation
cat("=== Creating Basic Manhattan Plot ===\n")
basic_plot <- create_manhattan_plot(
  con, 
  y_values = rda.simple.pq$q.values,
  vcf_file_id = 1,
  y_label = "RDA q-value",
  plot_title = "Basic RDA Manhattan Plot",
  verbose = TRUE
)

# Save basic plot
ggsave("basic_manhattan.pdf", basic_plot, width = 10, height = 6)
cat("Saved basic_manhattan.pdf\n")

# 2. Get functional summary from enrichment analysis
cat("\n=== Getting Functional Summary ===\n")
enrich_orf_summary <- summarize_functional_loci(
  con, 
  enrichment_results = enrich_orf$enrichment_results$BP, 
  candidate_file_id = 2,
  blast_param_id = 1,
  significance_threshold = 0.1
)

print(enrich_orf_summary$overview)

# 3. Create functional Manhattan plot with GO term labels (default)
cat("\n=== Creating Functional Manhattan Plot with GO Term Labels ===\n")
functional_plot_go <- create_functional_manhattan_plot(
  con, 
  y_values = rda.simple.pq$q.values,
  vcf_file_id = 1,
  functional_summary = enrich_orf_summary,
  y_label = "RDA q-value",
  plot_title = "RDA Analysis with GO Term Annotations",
  label_type = "go_term",  # Default
  verbose = TRUE
)

# Save functional plot with GO terms
ggsave("functional_manhattan_go.pdf", functional_plot_go, width = 12, height = 6)
cat("Saved functional_manhattan_go.pdf\n")

# 4. Create functional Manhattan plot with gene name labels
cat("\n=== Creating Functional Manhattan Plot with Gene Name Labels ===\n")
functional_plot_genes <- create_functional_manhattan_plot(
  con, 
  y_values = rda.simple.pq$q.values,
  vcf_file_id = 1,
  functional_summary = enrich_orf_summary,
  y_label = "RDA q-value",
  plot_title = "RDA Analysis with Gene Name Annotations",
  label_type = "gene_name",
  verbose = TRUE
)

# Save functional plot with gene names
ggsave("functional_manhattan_genes.pdf", functional_plot_genes, width = 12, height = 6)
cat("Saved functional_manhattan_genes.pdf\n")

# 5. Create functional Manhattan plot with UniProt accession labels
cat("\n=== Creating Functional Manhattan Plot with UniProt Accession Labels ===\n")
functional_plot_uniprot <- create_functional_manhattan_plot(
  con, 
  y_values = rda.simple.pq$q.values,
  vcf_file_id = 1,
  functional_summary = enrich_orf_summary,
  y_label = "RDA q-value",
  plot_title = "RDA Analysis with UniProt Accession Annotations",
  label_type = "uniprot_accession",
  verbose = TRUE
)

# Save functional plot with UniProt accessions
ggsave("functional_manhattan_uniprot.pdf", functional_plot_uniprot, width = 12, height = 6)
cat("Saved functional_manhattan_uniprot.pdf\n")

# 6. Create comparison plots for different BLAST runs
if (exists("enrich_raw_summary")) {
  cat("\n=== Creating Comparison with Raw Sequence Annotations ===\n")
  functional_plot_raw <- create_functional_manhattan_plot(
    con, 
    y_values = rda.simple.pq$q.values,
    vcf_file_id = 1,
    functional_summary = enrich_raw_summary,
    y_label = "RDA q-value",
    plot_title = "RDA Analysis with Raw Sequence Annotations",
    label_type = "go_term",
    verbose = TRUE
  )
  
  ggsave("functional_manhattan_raw.pdf", functional_plot_raw, width = 12, height = 6)
  cat("Saved functional_manhattan_raw.pdf\n")
}

# 7. Customized plot with different colors and threshold
cat("\n=== Creating Customized Manhattan Plot ===\n")
custom_plot <- create_functional_manhattan_plot(
  con, 
  y_values = rda.simple.pq$q.values,
  vcf_file_id = 1,
  functional_summary = enrich_orf_summary,
  y_label = "RDA q-value",
  plot_title = "Custom RDA Manhattan Plot",
  label_type = "go_term",
  signif_threshold = 0.05,  # Different threshold
  highlight_color = "#FF6B6B",  # Custom highlight color
  chr_colors = c("#2E8B57", "#4682B4"),  # Custom chromosome colors
  point_size = 1.5,  # Larger points
  verbose = TRUE
)

ggsave("custom_manhattan.pdf", custom_plot, width = 12, height = 6)
cat("Saved custom_manhattan.pdf\n")

# 8. Check chromosome consolidation status
cat("\n=== Chromosome Consolidation Status ===\n")
main_chroms <- DBI::dbGetQuery(con, "
  SELECT value FROM metadata WHERE key = 'main_chromosomes'
")

if (nrow(main_chroms) > 0) {
  defined_chroms <- jsonlite::fromJSON(main_chroms$value[1])
  cat("Main chromosomes defined:", paste(defined_chroms, collapse = ", "), "\n")
  cat("Scaffolds will be grouped as 'U' in plots\n")
} else {
  cat("No main chromosomes defined. Use define_chromosomes() to set up consolidation\n")
  cat("Example: define_chromosomes(con, c('LG1', 'LG2', 'LG3', ...))\n")
}

# Summary
cat("\n=== Summary ===\n")
cat("Generated Manhattan plots:\n")
cat("1. basic_manhattan.pdf - Basic plot without functional annotation\n")
cat("2. functional_manhattan_go.pdf - Functional plot with GO term labels\n")
cat("3. functional_manhattan_genes.pdf - Functional plot with gene name labels\n")
cat("4. functional_manhattan_uniprot.pdf - Functional plot with UniProt labels\n")
cat("5. custom_manhattan.pdf - Customized plot with different styling\n")

if (exists("enrich_raw_summary")) {
  cat("6. functional_manhattan_raw.pdf - Raw sequence annotation comparison\n")
}

cat("\nKey improvements:\n")
cat("- Automatic chromosome consolidation using database metadata\n")
cat("- Functional loci automatically highlighted and labeled\n")
cat("- Multiple label types: GO terms, gene names, UniProt accessions\n")
cat("- Proper chromosome ordering (LG1, LG2, ..., LG24, U)\n")
cat("- Clean x-axis labels without overlap\n")

close_funseq_db(con)