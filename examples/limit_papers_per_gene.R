#!/usr/bin/env Rscript

# Example: Limiting papers per gene to exactly 3
library(funseqR)
library(xml2)

# Connect to database
project_dir <- "/Users/brau0037/Library/CloudStorage/GoogleDrive-pygmyperch@gmail.com/My Drive/projects/snapperSG/south_coast/annotation/funseq_annotationMS"
db_path <- file.path(project_dir, "funseq_project.db")
con <- connect_funseq_db(db_path)

# Define enriched GO terms
enriched_go_ids <- c(
  "GO:0043204",  # C:perikaryon
  "GO:0005576",  # C:extracellular region
  "GO:0031463",  # C:Cul3-RING ubiquitin ligase complex
  "GO:0050839",  # F:cell adhesion molecule binding
  "GO:0030246"   # F:carbohydrate binding
)

cat("=== Literature Search with Strict Paper Limit ===\n")

# Run literature search
literature_results <- get_literature_for_enriched_genes(
  con = con,
  candidate_file_id = 3,
  background_file_id = 1,
  enriched_go_ids = enriched_go_ids,
  species = "danio rerio",
  max_papers_per_gene = 3,
  verbose = TRUE
)

cat("Before limiting:\n")
cat("Total papers:", nrow(literature_results$literature), "\n")

# Show papers per gene before limiting
if (nrow(literature_results$literature) > 0) {
  papers_per_gene <- table(literature_results$literature$uniprot_accession)
  cat("Papers per gene before limiting:\n")
  for (i in 1:length(papers_per_gene)) {
    cat("  ", names(papers_per_gene)[i], ": ", papers_per_gene[i], " papers\n")
  }
}

# Apply strict limit of 3 papers per gene
if (nrow(literature_results$literature) > 0) {
  
  cat("\nApplying strict limit of 3 papers per gene...\n")
  
  limited_list <- list()
  
  for (acc in unique(literature_results$literature$uniprot_accession)) {
    gene_papers <- literature_results$literature[literature_results$literature$uniprot_accession == acc, ]
    
    # Sort by year (newest first) then by relevance (could add more sorting criteria)
    if ("year" %in% names(gene_papers)) {
      gene_papers$year_numeric <- as.numeric(gene_papers$year)
      gene_papers <- gene_papers[order(-gene_papers$year_numeric, na.last = TRUE), ]
      gene_papers$year_numeric <- NULL
    }
    
    # Keep only first 3 papers
    if (nrow(gene_papers) > 3) {
      gene_papers <- gene_papers[1:3, ]
    }
    
    limited_list[[acc]] <- gene_papers
  }
  
  # Combine limited results
  literature_results$literature <- do.call(rbind, limited_list)
  rownames(literature_results$literature) <- NULL
  
  cat("After limiting:\n")
  cat("Total papers:", nrow(literature_results$literature), "\n")
  
  # Show papers per gene after limiting
  papers_per_gene_limited <- table(literature_results$literature$uniprot_accession)
  cat("Papers per gene after limiting:\n")
  for (i in 1:length(papers_per_gene_limited)) {
    cat("  ", names(papers_per_gene_limited)[i], ": ", papers_per_gene_limited[i], " papers\n")
  }
  
  # Generate references with limited papers
  refs <- generate_reference_list(literature_results$literature, format = "apa")
  writeLines(refs, "limited_literature_references.txt")
  
  cat("\nLimited references written to: limited_literature_references.txt\n")
  cat("Each gene now has maximum 3 papers (newest first)\n")
  
} else {
  cat("No literature found\n")
}

close_funseq_db(con)

cat("\n=== Why Some Genes Had >3 Papers Originally ===\n")
cat("The enhanced search tries multiple strategies per gene:\n")
cat("  1. Gene name search (up to 3 papers)\n")
cat("  2. UniProt accession search (up to 3 papers)\n") 
cat("  3. Functional term search (up to 3 papers)\n")
cat("  4. Custom search combinations (up to 3 papers)\n")
cat("  Total possible: up to 12+ papers per gene\n\n")

cat("The strict limiting approach keeps the 3 most recent papers per gene.\n")