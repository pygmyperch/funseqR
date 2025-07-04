# Test script to verify ORA p_adjusted fix
# Run this after the user re-runs their ORA analysis

library(funseqR)
library(RSQLite)

# Function to test p_adjusted storage
test_ora_pvalue_storage <- function(db_path) {
  
  cat("=== Testing ORA p_adjusted Storage Fix ===\n")
  
  if (!file.exists(db_path)) {
    cat("Database not found:", db_path, "\n")
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
  
  # Check ORA results
  ora_count <- dbGetQuery(con, "SELECT COUNT(*) as count FROM ora_results")$count
  cat("Total ORA results:", ora_count, "\n")
  
  if (ora_count == 0) {
    cat("No ORA results found. Please run ORA analysis first.\n")
    dbDisconnect(con)
    return(FALSE)
  }
  
  # Check p_adjusted values
  cat("\n=== Checking p_adjusted values ===\n")
  pvalue_check <- dbGetQuery(con, "
    SELECT 
      analysis_id,
      COUNT(*) as total_results,
      MIN(CAST(p_adjusted AS REAL)) as min_padj,
      MAX(CAST(p_adjusted AS REAL)) as max_padj,
      AVG(CAST(p_adjusted AS REAL)) as avg_padj,
      SUM(CASE WHEN CAST(p_adjusted AS REAL) <= 0.1 THEN 1 ELSE 0 END) as significant_count,
      SUM(CASE WHEN CAST(p_adjusted AS REAL) >= 1.0 THEN 1 ELSE 0 END) as suspicious_count
    FROM ora_results 
    GROUP BY analysis_id
    ORDER BY analysis_id
  ")
  
  print(pvalue_check)
  
  # Check for the fix
  all_good <- TRUE
  for (i in 1:nrow(pvalue_check)) {
    row <- pvalue_check[i, ]
    
    # Check if p_adjusted values are reasonable (between 0 and 1)
    if (row$min_padj < 0 || row$max_padj > 1) {
      cat("❌ Analysis ID", row$analysis_id, ": p_adjusted values outside [0,1] range\n")
      all_good <- FALSE
    } else if (row$suspicious_count > 0) {
      cat("⚠️  Analysis ID", row$analysis_id, ": Found", row$suspicious_count, "p_adjusted values >= 1.0 (suspicious)\n")
      all_good <- FALSE
    } else if (row$min_padj >= 1.0) {
      cat("❌ Analysis ID", row$analysis_id, ": All p_adjusted values >= 1.0 (still broken)\n")
      all_good <- FALSE
    } else {
      cat("✅ Analysis ID", row$analysis_id, ": p_adjusted values look correct (", 
          round(row$min_padj, 4), " - ", round(row$max_padj, 4), ")\n")
    }
  }
  
  # Show specific significant results
  if (any(pvalue_check$significant_count > 0)) {
    cat("\n=== Significant results (FDR < 0.1) ===\n")
    sig_results <- dbGetQuery(con, "
      SELECT analysis_id, term_id, substr(term_name, 1, 50) as term_name,
             p_value, p_adjusted
      FROM ora_results 
      WHERE CAST(p_adjusted AS REAL) <= 0.1
      ORDER BY CAST(p_adjusted AS REAL)
    ")
    print(sig_results)
  }
  
  dbDisconnect(con)
  
  if (all_good) {
    cat("\n🎉 SUCCESS: p_adjusted values are now stored correctly!\n")
    return(TRUE)
  } else {
    cat("\n💥 ISSUE: p_adjusted values still have problems. May need to re-run ORA.\n")
    return(FALSE)
  }
}

# Test with common database paths
possible_paths <- c(
  "funseq_project.db",
  "analysis.db"
)

for (db_path in possible_paths) {
  if (file.exists(db_path)) {
    cat("Testing database:", db_path, "\n")
    test_ora_pvalue_storage(db_path)
    break
  }
}