#!/usr/bin/env Rscript

# Clear and straightforward GO enrichment reporting for snapper analysis
library(funseqR)

# Connect to database
project_dir <- "/Users/brau0037/Library/CloudStorage/GoogleDrive-pygmyperch@gmail.com/My Drive/projects/snapperSG/south_coast/annotation/funseq_annotationMS"
db_path <- file.path(project_dir, "funseq_project.db")
con <- connect_funseq_db(db_path)

cat("=== Clear GO Enrichment Results for Snapper Analysis ===\n")

# Run GO enrichment analysis
cat("Running GO enrichment analysis...\n")
results <- run_go_enrichment_workflow(
  con = con,
  project_id = 1,
  candidate_vcf_file = file.path(project_dir, "SA448_855.vcf"),
  background_file_id = 1,
  ontologies = c("BP", "MF", "CC"),
  min_genes = 3,
  store_results = FALSE,
  create_plots = FALSE,
  verbose = FALSE
)

candidate_file_id <- results$summary$candidate_file_id
background_file_id <- results$summary$background_file_id

cat("Analysis complete.\n")
cat("Candidate genes:", results$summary$foreground_genes, "\n")
cat("Background genes:", results$summary$background_genes, "\n\n")

# Export 1: Summary of significantly enriched GO terms (unique terms only)
cat("1. Creating GO term enrichment summary...\n")

all_enriched_terms <- list()
ontology_summary <- data.frame(
  ontology = character(),
  total_tested = integer(),
  highly_significant = integer(),
  significant = integer(),
  trending = integer(),
  stringsAsFactors = FALSE
)

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

# Export ontology summary
write.csv(ontology_summary, "snapper_ontology_summary.csv", row.names = FALSE)

# Export enriched terms summary (unique terms only)
if (length(all_enriched_terms) > 0) {
  enriched_summary <- do.call(rbind, all_enriched_terms)
  enriched_summary <- enriched_summary[order(enriched_summary$p_adjusted), ]
  
  # Select key columns for clarity
  enriched_clean <- enriched_summary[, c(
    "ontology_name", "go_id", "go_term", "foreground_count", "background_count",
    "fold_enrichment", "p_value", "p_adjusted", "significance_level"
  )]
  
  write.csv(enriched_clean, "snapper_enriched_go_terms.csv", row.names = FALSE)
  
  cat("  Exported unique enriched GO terms:", nrow(enriched_clean), "\n")
} else {
  cat("  No significantly enriched GO terms found\n")
}

# Export 2: Gene-level summary for enriched terms
cat("2. Creating gene-level summary for enriched terms...\n")

if (length(all_enriched_terms) > 0) {
  
  # Get all significant GO IDs
  sig_go_ids <- unique(enriched_summary$go_id)
  
  # Query for genes associated with ANY significant GO term
  gene_query <- paste(
    "SELECT DISTINCT",
    "  c.chromosome,",
    "  c.position,", 
    "  c.ref,",
    "  c.alt,",
    "  a.uniprot_accession,",
    "  a.gene_names,",
    "  a.entry_name,",
    "  gt.go_id,",
    "  gt.go_term,",
    "  gt.go_category",
    "FROM vcf_data c",
    "JOIN vcf_data r ON (c.chromosome = r.chromosome AND c.position = r.position)",
    "JOIN flanking_sequences fs ON r.vcf_id = fs.vcf_id",
    "JOIN blast_results br ON fs.flanking_id = br.flanking_id",
    "JOIN annotations a ON br.blast_result_id = a.blast_result_id",
    "JOIN go_terms gt ON a.annotation_id = gt.annotation_id",
    "WHERE c.file_id = ? AND r.file_id = ?",
    "  AND gt.go_id IN (", paste0("'", sig_go_ids, "'", collapse = ", "), ")",
    "ORDER BY c.chromosome, c.position"
  )
  
  enriched_genes <- DBI::dbGetQuery(con, gene_query, 
                                   list(candidate_file_id, background_file_id))
  
  if (nrow(enriched_genes) > 0) {
    
    # Create gene-centric summary
    gene_keys <- paste(enriched_genes$chromosome, enriched_genes$position, 
                      enriched_genes$uniprot_accession, sep = "_")
    
    gene_summary_list <- list()
    for (key in unique(gene_keys)) {
      subset_data <- enriched_genes[gene_keys == key, ]
      
      # Get the enrichment stats for this gene's GO terms
      gene_go_ids <- unique(subset_data$go_id)
      gene_enrichment <- enriched_summary[enriched_summary$go_id %in% gene_go_ids, ]
      
      gene_summary_list[[key]] <- data.frame(
        chromosome = subset_data$chromosome[1],
        position = subset_data$position[1],
        ref_allele = subset_data$ref[1],
        alt_allele = subset_data$alt[1],
        uniprot_accession = subset_data$uniprot_accession[1],
        gene_names = subset_data$gene_names[1],
        entry_name = subset_data$entry_name[1],
        n_enriched_go_terms = length(gene_go_ids),
        best_p_adjusted = min(gene_enrichment$p_adjusted),
        max_fold_enrichment = max(gene_enrichment$fold_enrichment),
        enriched_go_terms = paste(unique(subset_data$go_term), collapse = " | "),
        stringsAsFactors = FALSE
      )
    }
    
    gene_summary <- do.call(rbind, gene_summary_list)
    gene_summary <- gene_summary[order(gene_summary$best_p_adjusted), ]
    
    write.csv(gene_summary, "snapper_genes_with_enriched_go.csv", row.names = FALSE)
    
    cat("  Exported genes with enriched GO terms:", nrow(gene_summary), "\n")
  }
}

# Export 3: Generate formatted text summary
cat("3. Creating formatted text summary...\n")

# Source the template functions
source("R/go_summary_template.R")

# Generate text summary
if (length(all_enriched_terms) > 0) {
  summary_text <- generate_go_summary_text(
    ontology_summary = ontology_summary,
    enriched_terms = enriched_summary,
    gene_count = if(exists("gene_summary")) nrow(gene_summary) else 0,
    output_file = "snapper_go_enrichment_summary.txt",
    format = "text"
  )
  
  # Also generate markdown version for reports
  markdown_summary <- generate_go_summary_text(
    ontology_summary = ontology_summary,
    enriched_terms = enriched_summary, 
    gene_count = if(exists("gene_summary")) nrow(gene_summary) else 0,
    output_file = "snapper_go_enrichment_summary.md",
    format = "markdown"
  )
  
  # Generate report content structure
  report_content <- generate_go_report_content(
    ontology_summary = ontology_summary,
    enriched_terms = enriched_summary,
    gene_count = if(exists("gene_summary")) nrow(gene_summary) else 0
  )
  
  # Save report content as RDS for dynamic report integration
  saveRDS(report_content, "snapper_go_report_content.rds")
  
  cat("  Generated formatted summaries and report content\n")
} else {
  # Handle case with no enriched terms
  summary_lines <- c(
    "GO Enrichment Analysis Summary",
    paste(rep("=", 35), collapse = ""),
    "",
    "Testing Overview:",
    ""
  )
  
  for (i in 1:nrow(ontology_summary)) {
    ont <- ontology_summary[i, ]
    ont_name <- switch(ont$ontology,
      "BP" = "Biological Process",
      "MF" = "Molecular Function", 
      "CC" = "Cellular Component"
    )
    
    line <- sprintf("  - %s: %d terms tested → %d enriched (FDR < 0.1)",
                   ont_name, ont$total_tested, ont$trending)
    summary_lines <- c(summary_lines, line)
  }
  
  summary_lines <- c(summary_lines, "", "No significantly enriched GO terms found (FDR < 0.1)")
  
  writeLines(summary_lines, "snapper_go_enrichment_summary.txt")
}

close_funseq_db(con)

cat("\n=== Clear Reporting Complete ===\n")
cat("Files generated:\n")

if (file.exists("snapper_ontology_summary.csv")) {
  cat("  snapper_ontology_summary.csv: Overview of testing results by ontology\n")
}

if (file.exists("snapper_enriched_go_terms.csv")) {
  cat("  snapper_enriched_go_terms.csv: Unique enriched GO terms (one row per term)\n")
}

if (file.exists("snapper_genes_with_enriched_go.csv")) {
  cat("  snapper_genes_with_enriched_go.csv: Genes contributing to enrichment\n")
}

if (file.exists("snapper_go_enrichment_summary.txt")) {
  cat("  snapper_go_enrichment_summary.txt: Human-readable text summary\n")
}

if (file.exists("snapper_go_enrichment_summary.md")) {
  cat("  snapper_go_enrichment_summary.md: Markdown summary for reports\n")
}

if (file.exists("snapper_go_report_content.rds")) {
  cat("  snapper_go_report_content.rds: Structured content for dynamic reports\n")
}

cat("\nInterpretation:\n")
cat("  - The ontology summary shows testing overview\n")
cat("  - The enriched terms file shows UNIQUE GO terms (avoid double-counting)\n") 
cat("  - The genes file shows which specific genes drive the enrichment\n")
cat("  - The text/markdown summaries provide human-readable interpretations\n")
cat("  - The RDS file contains structured data for dynamic report integration\n")