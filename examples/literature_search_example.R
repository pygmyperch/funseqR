#!/usr/bin/env Rscript

# Complete Literature Search Workflow Example
# This script demonstrates how to use the new literature search functionality

library(funseqR)
library(xml2)  # Required for PubMed API

# Step 1: Connect to database
project_dir <- "/Users/brau0037/Library/CloudStorage/GoogleDrive-pygmyperch@gmail.com/My Drive/projects/snapperSG/south_coast/annotation/funseq_annotationMS"
db_path <- file.path(project_dir, "funseq_project.db")
con <- connect_funseq_db(db_path)

cat("=== Literature Search Workflow Example ===\n")

# Step 2: Define enriched GO terms (from your GO enrichment analysis)
enriched_go_ids <- c(
  "GO:0043204",  # C:perikaryon
  "GO:0005576",  # C:extracellular region
  "GO:0031463",  # C:Cul3-RING ubiquitin ligase complex
  "GO:0050839",  # F:cell adhesion molecule binding
  "GO:0030246"   # F:carbohydrate binding
)

# Step 3: Run basic literature search
cat("Running basic literature search...\n")
literature_results <- get_literature_for_enriched_genes(
  con = con,
  candidate_file_id = 3,  # SA448_855.vcf
  background_file_id = 1, # SA448_14699.vcf
  enriched_go_ids = enriched_go_ids,
  species = "danio rerio",
  max_papers_per_gene = 3,
  verbose = TRUE
)

# Step 4: Display results summary
cat("\n=== Search Results Summary ===\n")
cat("Total genes searched:", literature_results$summary$total_genes, "\n")
cat("Papers found:", literature_results$summary$papers_found, "\n")
cat("Genes with literature:", literature_results$summary$genes_with_literature, "\n")

if (nrow(literature_results$literature) > 0) {
  
  # Step 5: Show which genes have literature
  cat("\n=== Genes with Literature Found ===\n")
  for (acc in unique(literature_results$literature$uniprot_accession)) {
    gene_papers <- literature_results$literature[literature_results$literature$uniprot_accession == acc, ]
    gene_info <- literature_results$genes[literature_results$genes$uniprot_accession == acc, ]
    
    gene_display <- if(gene_info$gene_names[1] != "") {
      paste0(acc, " (", gene_info$gene_names[1], ")")
    } else {
      paste0(acc, " (", gene_info$entry_name[1], ")")
    }
    
    cat("  ", gene_display, ": ", nrow(gene_papers), " papers\n")
    
    # Show first paper as example
    if (nrow(gene_papers) > 0) {
      cat("    Example: ", gene_papers$title[1], " (", gene_papers$year[1], ")\n")
    }
  }
  
  # Step 6: Generate reference files
  cat("\n=== Generating Reference Files ===\n")
  
  # APA format references
  apa_refs <- generate_reference_list(
    literature_results$literature,
    output_file = "my_enriched_genes_references_apa.txt",
    format = "apa"
  )
  
  # Detailed format with abstracts and search info
  detailed_refs <- generate_reference_list(
    literature_results$literature,
    output_file = "my_enriched_genes_references_detailed.txt",
    format = "detailed"
  )
  
  # Export complete data
  write.csv(literature_results$literature, "my_literature_search_complete.csv", row.names = FALSE)
  write.csv(literature_results$genes, "my_enriched_genes_info.csv", row.names = FALSE)
  
  cat("Files generated:\n")
  cat("  my_enriched_genes_references_apa.txt\n")
  cat("  my_enriched_genes_references_detailed.txt\n")
  cat("  my_literature_search_complete.csv\n")
  cat("  my_enriched_genes_info.csv\n")
  
  # Step 7: Show search strategy effectiveness
  cat("\n=== Most Effective Search Terms ===\n")
  search_summary <- table(literature_results$literature$search_term)
  for (term in names(search_summary)) {
    cat("  '", term, "': ", search_summary[term], " papers\n", sep = "")
  }
  
} else {
  cat("\nNo literature found. Try:\n")
  cat("  - Different custom search terms\n")
  cat("  - Broader species terms\n") 
  cat("  - Manual searches for key genes\n")
}

close_funseq_db(con)

cat("\n=== Example: Enhanced Search with Custom Terms ===\n")
cat("# For snapper-specific research:\n")
cat("literature <- get_literature_for_enriched_genes(\n")
cat("  con, 3, 1, enriched_go_ids,\n")
cat("  species = 'chrysophrys auratus',\n")
cat("  custom_search_term = 'snapper'\n")
cat(")\n\n")

cat("# For broader marine research:\n")
cat("literature <- get_literature_for_enriched_genes(\n")
cat("  con, 3, 1, enriched_go_ids,\n")
cat("  species = 'teleost',\n")
cat("  custom_search_term = 'marine adaptation'\n")
cat(")\n\n")

cat("=== Literature Search Workflow Complete ===\n")