#!/usr/bin/env Rscript

# Debug the exact GO enrichment storage issue
library(funseqR)

project_dir <- "/Users/brau0037/Library/CloudStorage/GoogleDrive-pygmyperch@gmail.com/My Drive/projects/snapperSG/south_coast/annotation/funseq_annotationMS"
db_path <- file.path(project_dir, "funseq_project.db")

con <- connect_funseq_db(db_path)

# Create mock results similar to your structure
results <- data.frame(
  go_id = c("GO:0031640", "GO:0006397", "GO:0008380"),
  go_term = c("P:killing of cells of another organism", "P:mRNA processing", "P:RNA splicing"),
  go_category = c("BP", "BP", "BP"),
  foreground_count = c(5, 4, 3),
  background_count = c(9, 7, 5),
  total_foreground = c(124, 124, 124),
  total_background = c(1176, 1176, 1176),
  expected_count = c(0.949, 0.738, 0.527),
  fold_enrichment = c(5.27, 5.42, 5.69),
  p_value = c(0.00107, 0.00321, 0.00977),
  p_adjusted = c(0.131, 0.196, 0.254),
  significance_level = c("not_significant", "not_significant", "not_significant"),
  stringsAsFactors = FALSE
)

project_id <- 1
candidate_file_id <- 3
background_file_id <- 1
ontology <- "BP"
parameters <- list(min_genes = 3, max_genes = 500)
verbose <- TRUE

cat("Testing store_go_enrichment_results function...\n")

tryCatch({
  enrichment_id <- store_go_enrichment_results(
    con, 
    project_id, 
    candidate_file_id, 
    background_file_id, 
    results, 
    ontology, 
    parameters = parameters, 
    verbose = verbose
  )
  cat("SUCCESS: Function completed with enrichment_id:", enrichment_id, "\n")
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  cat("Full error object:\n")
  print(e)
})

close_funseq_db(con)