#!/usr/bin/env Rscript

# Debug step by step
library(funseqR)

# Replicate the exact function logic step by step
project_dir <- "/Users/brau0037/Library/CloudStorage/GoogleDrive-pygmyperch@gmail.com/My Drive/projects/snapperSG/south_coast/annotation/funseq_annotationMS"
db_path <- file.path(project_dir, "funseq_project.db")

con <- connect_funseq_db(db_path)

# Mock data
enrichment_results <- data.frame(
  go_id = c("GO:0031640", "GO:0006397"),
  go_term = c("P:killing of cells", "P:mRNA processing"),
  go_category = c("BP", "BP"),
  foreground_count = c(5, 4),
  background_count = c(9, 7),
  total_foreground = c(124, 124),
  total_background = c(1176, 1176),
  expected_count = c(0.949, 0.738),
  fold_enrichment = c(5.27, 5.42),
  p_value = c(0.00107, 0.00321),
  p_adjusted = c(0.131, 0.196),
  significance_level = c("not_significant", "not_significant"),
  stringsAsFactors = FALSE
)

project_id <- 1
foreground_file_id <- 3
background_file_id <- 1
ontology <- "BP"
parameters <- list(min_genes = 3, max_genes = 500)

cat("Step 1: Checking tables...\n")
tables <- DBI::dbListTables(con)
if (!all(c("go_enrichment_analyses", "go_enrichment_results") %in% tables)) {
  stop("GO enrichment tables not found in database. Please upgrade schema.")
}
cat("Tables exist: OK\n")

cat("Step 2: Preparing parameters...\n")
if (is.null(parameters)) {
  parameters <- list()
}
params_json <- jsonlite::toJSON(parameters, auto_unbox = TRUE)
cat("Parameters JSON created: OK\n")

cat("Step 3: Testing analysis insert...\n")
analysis_query <- "
  INSERT INTO go_enrichment_analyses 
  (project_id, foreground_file_id, background_file_id, ontology, analysis_date, 
   total_foreground_genes, total_background_genes, analysis_parameters)
  VALUES (?, ?, ?, ?, ?, ?, ?, ?)
"

tryCatch({
  DBI::dbExecute(con, analysis_query, list(
    project_id,
    foreground_file_id,
    background_file_id,
    ontology,
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    if(nrow(enrichment_results) > 0) enrichment_results$total_foreground[1] else 0,
    if(nrow(enrichment_results) > 0) enrichment_results$total_background[1] else 0,
    params_json
  ))
  cat("Analysis insert: SUCCESS\n")
}, error = function(e) {
  cat("Analysis insert: ERROR -", e$message, "\n")
  stop("Analysis insert failed")
})

# Get the analysis ID
enrichment_id <- DBI::dbGetQuery(con, "SELECT last_insert_rowid() as id")$id
cat("Enrichment ID:", enrichment_id, "\n")

cat("Step 4: Testing results insert...\n")
if (nrow(enrichment_results) > 0) {
  
  # Prepare results for insertion
  results_to_insert <- enrichment_results
  results_to_insert$enrichment_id <- enrichment_id
  
  # Select only the columns we need in the right order
  results_to_insert <- results_to_insert[, c("enrichment_id", "go_id", "go_term", "go_category", 
                                           "foreground_count", "background_count", "total_foreground", 
                                           "total_background", "expected_count", "fold_enrichment", 
                                           "p_value", "p_adjusted", "significance_level")]
  
  cat("Prepared data structure:\n")
  str(results_to_insert)
  
  result_query <- "
    INSERT INTO go_enrichment_results 
    (enrichment_id, go_id, go_term, go_category, foreground_count, background_count,
     total_foreground, total_background, expected_count, fold_enrichment, 
     p_value, p_adjusted, significance_level)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
  "
  
  cat("Testing first row insert...\n")
  row_data <- as.list(results_to_insert[1, ])
  cat("Row data structure:\n")
  str(row_data)
  
  tryCatch({
    DBI::dbExecute(con, result_query, row_data)
    cat("First row insert: SUCCESS\n")
  }, error = function(e) {
    cat("First row insert: ERROR -", e$message, "\n")
  })
}

close_funseq_db(con)