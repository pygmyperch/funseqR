#!/usr/bin/env Rscript

# Test the exact database insert that's failing
library(funseqR)

project_dir <- "/Users/brau0037/Library/CloudStorage/GoogleDrive-pygmyperch@gmail.com/My Drive/projects/snapperSG/south_coast/annotation/funseq_annotationMS"
db_path <- file.path(project_dir, "funseq_project.db")

con <- connect_funseq_db(db_path)

cat("Testing GO enrichment analysis insert...\n")

analysis_query <- "
  INSERT INTO go_enrichment_analyses 
  (project_id, foreground_file_id, background_file_id, ontology, analysis_date, 
   total_foreground_genes, total_background_genes, analysis_parameters)
  VALUES (?, ?, ?, ?, ?, ?, ?, ?)
"

params_json <- jsonlite::toJSON(list(min_genes = 3, max_genes = 500), auto_unbox = TRUE)

tryCatch({
  DBI::dbExecute(con, analysis_query, list(
    1,  # project_id
    3,  # foreground_file_id  
    1,  # background_file_id
    "TEST",  # ontology
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),  # analysis_date
    124,  # total_foreground_genes
    1176,  # total_background_genes
    params_json  # analysis_parameters
  ))
  cat("SUCCESS: Analysis insert worked!\n")
}, error = function(e) {
  cat("ERROR in analysis insert:", e$message, "\n")
})

close_funseq_db(con)