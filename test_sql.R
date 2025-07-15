#!/usr/bin/env Rscript

# Test the exact SQL query that's failing
library(funseqR)

cat("=== Testing SQL Query Directly ===\n")

# Connect to database
con <- connect_funseq_db('funseq_project.db')

# Test the exact same query with simple values
analysis_query <- "
  INSERT INTO ora_analyses
  (blast_param_id, annotation_type, term_type, analysis_date,
   total_foreground_genes, total_background_genes, analysis_parameters, enrichment_method)
  VALUES (?, ?, ?, ?, ?, ?, ?, ?)
"

# Test parameters (exactly as they would be in the function)
# For "one database = one analysis" model, always use blast_param_id = 1
db_blast_param_id <- 1L
db_annotation_type <- "GO"
db_term_type <- "BP"
timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
total_fg <- 39L
total_bg <- 506L
params_json <- "{}"
db_method <- "clusterprofiler"

cat("Testing SQL execution with these parameters:\n")
cat("  1. db_blast_param_id:", if(is.null(db_blast_param_id)) "NULL" else db_blast_param_id, "\n")
cat("  2. db_annotation_type:", db_annotation_type, "\n")
cat("  3. db_term_type:", db_term_type, "\n")
cat("  4. timestamp:", timestamp, "\n")
cat("  5. total_fg:", total_fg, "\n")
cat("  6. total_bg:", total_bg, "\n")
cat("  7. params_json:", params_json, "\n")
cat("  8. db_method:", db_method, "\n")

tryCatch({
  result <- DBI::dbExecute(con, analysis_query, list(
    db_blast_param_id,
    db_annotation_type,
    db_term_type,
    timestamp,
    total_fg,
    total_bg,
    params_json,
    db_method
  ))
  cat("SUCCESS: SQL executed, rows affected:", result, "\n")
  
  # Get the analysis ID
  analysis_id <- DBI::dbGetQuery(con, "SELECT last_insert_rowid() as id")$id
  cat("Analysis ID created:", analysis_id, "\n")
  
}, error = function(e) {
  cat("ERROR in SQL execution:", e$message, "\n")
  cat("Error class:", class(e), "\n")
  print(e)
})

# Close connection
close_funseq_db(con)