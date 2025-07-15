#!/usr/bin/env Rscript

# Update schema to support stored candidates workflow
library(funseqR)

update_ora_schema <- function(con, verbose = TRUE) {
  if (verbose) message("Updating ORA schema to support stored candidates...")
  
  # Check if ora_analyses table exists
  tables <- DBI::dbListTables(con)
  if (!"ora_analyses" %in% tables) {
    if (verbose) message("ora_analyses table doesn't exist, nothing to update")
    return(invisible(NULL))
  }
  
  # Get current schema
  schema_info <- DBI::dbGetQuery(con, "PRAGMA table_info(ora_analyses)")
  
  # Check if foreground_file_id is currently NOT NULL
  fg_col <- schema_info[schema_info$name == "foreground_file_id", ]
  if (nrow(fg_col) > 0 && fg_col$notnull == 1) {
    if (verbose) message("Updating foreground_file_id to allow NULL values...")
    
    # Begin transaction
    DBI::dbExecute(con, "BEGIN TRANSACTION")
    
    tryCatch({
      # Create new table with correct schema
      DBI::dbExecute(con, "
        CREATE TABLE ora_analyses_new (
          analysis_id INTEGER PRIMARY KEY,
          foreground_file_id INTEGER,
          background_file_id INTEGER NOT NULL,
          blast_param_id INTEGER,
          annotation_type TEXT NOT NULL,
          term_type TEXT NOT NULL,
          analysis_date TEXT NOT NULL,
          total_foreground_genes INTEGER,
          total_background_genes INTEGER,
          analysis_parameters TEXT,
          enrichment_method TEXT DEFAULT 'clusterprofiler',
          FOREIGN KEY (foreground_file_id) REFERENCES input_files (file_id),
          FOREIGN KEY (background_file_id) REFERENCES input_files (file_id),
          FOREIGN KEY (blast_param_id) REFERENCES blast_parameters (blast_param_id)
        )
      ")
      
      # Copy existing data
      existing_data <- DBI::dbGetQuery(con, "SELECT COUNT(*) as count FROM ora_analyses")
      if (existing_data$count > 0) {
        DBI::dbExecute(con, "
          INSERT INTO ora_analyses_new 
          SELECT * FROM ora_analyses
        ")
        if (verbose) message("Copied ", existing_data$count, " existing records")
      }
      
      # Drop old table and rename new one
      DBI::dbExecute(con, "DROP TABLE ora_analyses")
      DBI::dbExecute(con, "ALTER TABLE ora_analyses_new RENAME TO ora_analyses")
      
      # Recreate indexes
      DBI::dbExecute(con, "CREATE INDEX idx_ora_fg_file ON ora_analyses (foreground_file_id)")
      DBI::dbExecute(con, "CREATE INDEX idx_ora_bg_file ON ora_analyses (background_file_id)")
      DBI::dbExecute(con, "CREATE INDEX idx_ora_blast_param ON ora_analyses (blast_param_id)")
      DBI::dbExecute(con, "CREATE INDEX idx_ora_annotation_type ON ora_analyses (annotation_type)")
      DBI::dbExecute(con, "CREATE INDEX idx_ora_term_type ON ora_analyses (term_type)")
      DBI::dbExecute(con, "CREATE INDEX idx_ora_method ON ora_analyses (enrichment_method)")
      
      # Commit transaction
      DBI::dbExecute(con, "COMMIT")
      
      if (verbose) message("Schema update completed successfully")
      
    }, error = function(e) {
      # Rollback on error
      DBI::dbExecute(con, "ROLLBACK")
      stop("Schema update failed: ", e$message)
    })
    
  } else {
    if (verbose) message("Schema already supports NULL foreground_file_id")
  }
  
  return(invisible(NULL))
}

# Apply the update
con <- connect_funseq_db('funseq_project.db')
update_ora_schema(con, verbose = TRUE)
close_funseq_db(con)

cat("Schema update completed!\n")