#!/usr/bin/env Rscript

# Fixed Literature Search Example
library(funseqR)
library(xml2)

# Connect to database
project_dir <- "/Users/brau0037/Library/CloudStorage/GoogleDrive-pygmyperch@gmail.com/My Drive/projects/snapperSG/south_coast/annotation/funseq_annotationMS"
db_path <- file.path(project_dir, "funseq_project.db")
con <- connect_funseq_db(db_path)

cat("=== Fixed Literature Search Example ===\n")

# Define enriched GO terms
enriched_go_ids <- c(
  "GO:0043204",  # C:perikaryon
  "GO:0005576",  # C:extracellular region
  "GO:0031463",  # C:Cul3-RING ubiquitin ligase complex
  "GO:0050839",  # F:cell adhesion molecule binding
  "GO:0030246"   # F:carbohydrate binding
)

# Run literature search with custom search term
cat("Running enhanced literature search...\n")
literature_results <- get_literature_for_enriched_genes(
  con = con,
  candidate_file_id = 3,
  background_file_id = 1,
  enriched_go_ids = enriched_go_ids,
  species = "danio rerio",
  custom_search_term = "snapper",  # Custom search term
  max_papers_per_gene = 3,
  verbose = TRUE
)

cat("\n=== Results Summary ===\n")
cat("Total genes:", literature_results$summary$total_genes, "\n")
cat("Papers found:", literature_results$summary$papers_found, "\n")
cat("Genes with literature:", literature_results$summary$genes_with_literature, "\n")

if (nrow(literature_results$literature) > 0) {
  
  # Generate fixed reference list
  cat("\nGenerating fixed reference list...\n")
  refs <- generate_reference_list(
    literature_results$literature,
    format = "apa"
  )
  
  # Write to file
  writeLines(refs, "fixed_literature_references.txt")
  cat("References written to: fixed_literature_references.txt\n")
  
  # Also save raw data
  write.csv(literature_results$literature, "fixed_literature_data.csv", row.names = FALSE)
  cat("Raw data saved to: fixed_literature_data.csv\n")
  
  # Show structure of results
  cat("\nData structure:\n")
  cat("Columns:", paste(names(literature_results$literature), collapse = ", "), "\n")
  cat("Number of references:", length(refs), "\n")
  
  # Show first few lines of references
  cat("\nFirst few reference lines:\n")
  cat(paste(head(refs, 10), collapse = ""))
  
} else {
  cat("No literature found\n")
}

close_funseq_db(con)

cat("\n=== How to Write References to File ===\n")
cat("# Method 1: writeLines() - recommended\n")
cat("writeLines(refs, 'my_references.txt')\n\n")

cat("# Method 2: cat() with formatting\n")
cat("cat(refs, file = 'my_references.txt', sep = '\\n')\n\n")

cat("# Method 3: with timestamp\n")
cat("timestamp <- format(Sys.time(), '%Y%m%d_%H%M')\n")
cat("filename <- paste0('references_', timestamp, '.txt')\n")
cat("writeLines(refs, filename)\n\n")

cat("=== Fixed Literature Search Complete ===\n")