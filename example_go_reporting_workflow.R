#!/usr/bin/env Rscript

# Complete GO Enrichment Analysis and Reporting Workflow
# This script demonstrates the full pipeline from analysis to formatted reporting

library(funseqR)

# Database connection
project_dir <- "/Users/brau0037/Library/CloudStorage/GoogleDrive-pygmyperch@gmail.com/My Drive/projects/snapperSG/south_coast/annotation/funseq_annotationMS"
db_path <- file.path(project_dir, "funseq_project.db")
con <- connect_funseq_db(db_path)

cat("=== Complete GO Enrichment Analysis and Reporting Workflow ===\n")

# Step 1: Run GO enrichment analysis
cat("Step 1: Running GO enrichment analysis...\n")
results <- run_go_enrichment_workflow(
  con = con,
  project_id = 1,
  candidate_vcf_file = file.path(project_dir, "SA448_855.vcf"),
  background_file_id = 1,
  ontologies = c("BP", "MF", "CC"),
  min_genes = 3,
  store_results = TRUE,  # Store in database
  create_plots = FALSE,
  verbose = FALSE
)

cat("  Analysis complete: ", results$summary$foreground_genes, " candidate vs ", 
    results$summary$background_genes, " background genes\n")

# Step 2: Generate clear summary files
cat("Step 2: Generating clear summary files...\n")

# Create ontology summary
ontology_summary <- data.frame(
  ontology = character(),
  total_tested = integer(),
  highly_significant = integer(),
  significant = integer(),
  trending = integer(),
  stringsAsFactors = FALSE
)

all_enriched_terms <- list()
for (ontology in c("BP", "MF", "CC")) {
  if (!is.null(results$enrichment_results[[ontology]])) {
    
    go_results <- results$enrichment_results[[ontology]]
    
    # Count significance levels
    highly_sig <- sum(go_results$p_adjusted < 0.01)
    sig <- sum(go_results$p_adjusted < 0.05) 
    trending <- sum(go_results$p_adjusted < 0.10)
    
    ontology_summary <- rbind(ontology_summary, data.frame(
      ontology = ontology,
      total_tested = nrow(go_results),
      highly_significant = highly_sig,
      significant = sig,
      trending = trending,
      stringsAsFactors = FALSE
    ))
    
    # Get significant/trending terms
    sig_terms <- go_results[go_results$p_adjusted < 0.10, ]
    
    if (nrow(sig_terms) > 0) {
      sig_terms$ontology_name <- switch(ontology,
        "BP" = "Biological Process",
        "MF" = "Molecular Function", 
        "CC" = "Cellular Component"
      )
      all_enriched_terms[[ontology]] <- sig_terms
    }
  }
}

# Combine enriched terms
if (length(all_enriched_terms) > 0) {
  enriched_summary <- do.call(rbind, all_enriched_terms)
  enriched_summary <- enriched_summary[order(enriched_summary$p_adjusted), ]
} else {
  enriched_summary <- data.frame()
}

# Step 3: Generate gene-level summary
cat("Step 3: Creating gene-level summary...\n")
gene_count <- 0

if (nrow(enriched_summary) > 0) {
  # Get genes associated with enriched terms
  sig_go_ids <- unique(enriched_summary$go_id)
  
  gene_query <- paste(
    "SELECT DISTINCT",
    "  c.chromosome,",
    "  c.position,", 
    "  c.ref,",
    "  c.alt,",
    "  a.uniprot_accession,",
    "  a.gene_names",
    "FROM vcf_data c",
    "JOIN vcf_data r ON (c.chromosome = r.chromosome AND c.position = r.position)",
    "JOIN flanking_sequences fs ON r.vcf_id = fs.vcf_id",
    "JOIN blast_results br ON fs.flanking_id = br.flanking_id",
    "JOIN annotations a ON br.blast_result_id = a.blast_result_id",
    "JOIN go_terms gt ON a.annotation_id = gt.annotation_id",
    "WHERE c.file_id = ? AND r.file_id = ?",
    "  AND gt.go_id IN (", paste0("'", sig_go_ids, "'", collapse = ", "), ")"
  )
  
  enriched_genes <- DBI::dbGetQuery(con, gene_query, 
                                   list(results$summary$candidate_file_id, 
                                        results$summary$background_file_id))
  
  if (nrow(enriched_genes) > 0) {
    gene_keys <- paste(enriched_genes$chromosome, enriched_genes$position, 
                      enriched_genes$uniprot_accession, sep = "_")
    gene_count <- length(unique(gene_keys))
  }
}

# Step 4: Generate formatted summaries using templates
cat("Step 4: Generating formatted summaries...\n")

# Source template functions
source("R/go_summary_template.R")

# Generate text summary
text_summary <- generate_go_summary_text(
  ontology_summary = ontology_summary,
  enriched_terms = enriched_summary,
  gene_count = gene_count,
  output_file = "workflow_go_summary.txt",
  format = "text"
)

# Generate markdown summary
markdown_summary <- generate_go_summary_text(
  ontology_summary = ontology_summary,
  enriched_terms = enriched_summary, 
  gene_count = gene_count,
  output_file = "workflow_go_summary.md",
  format = "markdown"
)

# Generate structured report content
report_content <- generate_go_report_content(
  ontology_summary = ontology_summary,
  enriched_terms = enriched_summary,
  gene_count = gene_count
)

# Save structured content
saveRDS(report_content, "workflow_go_report_content.rds")

# Step 5: Add to dynamic report (if report exists)
cat("Step 5: Integrating with dynamic report...\n")

# Check if project has a report
report_info <- get_report_info(con, 1)

if (!is.null(report_info)) {
  
  cat("  Found existing report, adding GO enrichment section...\n")
  
  # Add GO enrichment to the report
  updated_report <- add_go_enrichment_to_report(
    con = con,
    project_id = 1,
    go_content_file = "workflow_go_report_content.rds",
    verbose = TRUE
  )
  
} else {
  cat("  No existing report found. Creating new report...\n")
  
  # Create a new report first
  report_file <- create_analysis_report(
    con = con,
    project_id = 1,
    output_file = "snapper_analysis_report.Rmd",
    template = "standard",
    render = FALSE
  )
  
  # Then add GO enrichment section
  updated_report <- add_go_enrichment_to_report(
    con = con,
    project_id = 1,
    go_content_file = "workflow_go_report_content.rds",
    verbose = TRUE
  )
}

close_funseq_db(con)

cat("\n=== Complete Workflow Summary ===\n")
cat("Generated files:\n")
files_to_check <- c(
  "workflow_go_summary.txt",
  "workflow_go_summary.md", 
  "workflow_go_report_content.rds"
)

for (file in files_to_check) {
  if (file.exists(file)) {
    cat("  ✓", file, "\n")
  } else {
    cat("  ✗", file, "(not found)\n")
  }
}

cat("\nDynamic report:\n")
if (exists("updated_report")) {
  cat("  ✓ Report updated:", updated_report, "\n")
} else {
  cat("  ✗ Report integration failed\n")
}

cat("\nWorkflow complete! The GO enrichment analysis has been:\n")
cat("  1. Performed and stored in the database\n")
cat("  2. Exported as clear, structured summaries\n") 
cat("  3. Formatted using templates for consistent reporting\n")
cat("  4. Integrated into the dynamic project report\n")

# Display summary
cat("\n" , paste(rep("=", 50), collapse = ""), "\n")
cat("FINAL SUMMARY\n")
cat(paste(rep("=", 50), collapse = ""), "\n")
cat(paste(text_summary, collapse = "\n"))
cat("\n", paste(rep("=", 50), collapse = ""), "\n")