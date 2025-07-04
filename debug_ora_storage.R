# Debug script to verify ORA p_adjusted storage fix
library(RSQLite)
library(DBI)

debug_ora_storage <- function(db_path = "funseq_project.db") {
  
  cat("=== ORA Storage Debug Script ===\n")
  
  if (!file.exists(db_path)) {
    cat("Database not found:\", db_path, \"\n")
    return(FALSE)
  }
  
  con <- dbConnect(SQLite(), db_path)
  
  # Check if ORA tables exist
  tables <- dbListTables(con)
  if (!all(c("ora_analyses", "ora_results") %in% tables)) {
    cat("ORA tables not found in database\n")
    dbDisconnect(con)
    return(FALSE)
  }
  
  # Check stored values
  cat("Checking stored p_adjusted values...\n")
  stored_check <- dbGetQuery(con, "
    SELECT analysis_id, term_id, 
           p_adjusted, 
           typeof(p_adjusted) as type,
           CASE WHEN p_adjusted LIKE '%/%' THEN 'FRACTION' ELSE 'DECIMAL' END as format
    FROM ora_results 
    ORDER BY analysis_id, CAST(p_adjusted AS REAL)
    LIMIT 10
  ")
  
  if (nrow(stored_check) > 0) {
    print(stored_check)
    
    # Check for significant terms
    cat("\nLooking for significant terms (p_adjusted < 0.1)...\n")
    significant_terms <- dbGetQuery(con, "
      SELECT analysis_id, term_id, term_name, p_adjusted,
             CAST(p_adjusted AS REAL) as numeric_padj
      FROM ora_results 
      WHERE CAST(p_adjusted AS REAL) < 0.1
      ORDER BY CAST(p_adjusted AS REAL)
    ")
    
    if (nrow(significant_terms) > 0) {
      cat("Found", nrow(significant_terms), "significant terms:\n")
      print(significant_terms)
    } else {
      cat("No significant terms found\n")
    }
  } else {
    cat("No ORA results found in database\n")
  }
  
  dbDisconnect(con)
  cat("=== Debug Complete ===\n")
  return(TRUE)
}

# Run the debug
debug_ora_storage()