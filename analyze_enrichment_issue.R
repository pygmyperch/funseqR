#!/usr/bin/env Rscript

# Targeted analysis of the enrichment extraction issue
# Based on the code analysis, here are the potential problems:

cat("=== Analysis of Potential Issues ===\n")

cat("\nPOTENTIAL ISSUE #1: Data Type Mismatch\n")
cat("The p_adjusted column might be stored as TEXT instead of REAL\n")
cat("Solution: Use CAST(p_adjusted AS REAL) in the comparison\n")

cat("\nPOTENTIAL ISSUE #2: NULL or Empty Values\n") 
cat("The p_adjusted column might contain NULL or empty string values\n")
cat("Solution: Add NULL checks in the WHERE clause\n")

cat("\nPOTENTIAL ISSUE #3: Floating Point Precision\n")
cat("Floating point comparisons might have precision issues\n") 
cat("Solution: Use a small epsilon value or ROUND function\n")

cat("\nPOTENTIAL ISSUE #4: Analysis ID Mismatch\n")
cat("The analysis_ids might not exist in the database\n")
cat("Solution: Verify the analysis_ids exist before querying\n")

cat("\nPOTENTIAL ISSUE #5: JOIN Issues\n")
cat("The JOIN between ora_analyses and ora_results might fail\n")
cat("Solution: Check foreign key relationships\n")

cat("\n=== RECOMMENDED FIXES ===\n")

cat("\nFIX #1: Improved Query with Data Type Handling\n")
improved_query <- "
SELECT DISTINCT
  ora.analysis_id,
  ora.annotation_type,
  ora.term_type,
  res.term_id,
  res.term_name,
  res.p_value,
  res.p_adjusted as fdr,
  res.gene_ids,
  res.fold_enrichment
FROM ora_analyses ora
JOIN ora_results res ON ora.analysis_id = res.analysis_id
WHERE ora.analysis_id IN (1, 2, 3)
  AND res.p_adjusted IS NOT NULL
  AND res.p_adjusted != ''
  AND CAST(res.p_adjusted AS REAL) <= 0.1
ORDER BY CAST(res.p_adjusted AS REAL)
"
cat(improved_query)

cat("\nFIX #2: Updated .extract_enrichment_data function\n")
cat("Replace the existing function with better error handling:\n")

improved_function <- '
.extract_enrichment_data <- function(con, analysis_ids, significance_threshold, verbose) {
  
  if (verbose) message("    - Retrieving enrichment results from database...")
  
  # First check if analysis_ids exist
  analysis_ids_str <- paste(analysis_ids, collapse = ", ")
  
  check_query <- paste0("
    SELECT analysis_id FROM ora_analyses 
    WHERE analysis_id IN (", analysis_ids_str, ")
  ")
  
  existing_ids <- DBI::dbGetQuery(con, check_query)$analysis_id
  
  if (length(existing_ids) == 0) {
    if (verbose) message("    - No analyses found for IDs: ", analysis_ids_str)
    return(data.frame())
  }
  
  if (length(existing_ids) < length(analysis_ids)) {
    missing_ids <- setdiff(analysis_ids, existing_ids)
    if (verbose) message("    - Warning: Analysis IDs not found: ", paste(missing_ids, collapse = ", "))
  }
  
  # Use only existing IDs
  analysis_ids_str <- paste(existing_ids, collapse = ", ")
  
  # Improved query with data type handling
  enrichment_query <- paste0("
    SELECT DISTINCT
      ora.analysis_id,
      ora.annotation_type,
      ora.term_type,
      res.term_id,
      res.term_name,
      res.p_value,
      res.p_adjusted as fdr,
      res.gene_ids,
      res.fold_enrichment
    FROM ora_analyses ora
    JOIN ora_results res ON ora.analysis_id = res.analysis_id
    WHERE ora.analysis_id IN (", analysis_ids_str, ")
      AND res.p_adjusted IS NOT NULL
      AND res.p_adjusted != \\\"\\\"
      AND CAST(res.p_adjusted AS REAL) <= ?
    ORDER BY CAST(res.p_adjusted AS REAL)
  ")
  
  enrichment_data <- DBI::dbGetQuery(con, enrichment_query, list(significance_threshold))
  
  if (verbose) {
    message("    - Found ", nrow(enrichment_data), " significantly enriched terms")
    if (nrow(enrichment_data) > 0) {
      by_analysis <- table(enrichment_data$analysis_id)
      for (i in names(by_analysis)) {
        message("      - Analysis ID ", i, ": ", by_analysis[i], " terms")
      }
    } else {
      # Debug output
      debug_query <- paste0("
        SELECT 
          ora.analysis_id,
          COUNT(*) as total_results,
          MIN(CAST(res.p_adjusted AS REAL)) as min_padj,
          MAX(CAST(res.p_adjusted AS REAL)) as max_padj
        FROM ora_analyses ora
        JOIN ora_results res ON ora.analysis_id = res.analysis_id
        WHERE ora.analysis_id IN (", analysis_ids_str, ")
          AND res.p_adjusted IS NOT NULL
        GROUP BY ora.analysis_id
      ")
      debug_info <- DBI::dbGetQuery(con, debug_query)
      message("    - Debug: Results by analysis ID:")
      print(debug_info)
      message("    - Significance threshold: ", significance_threshold)
    }
  }
  
  return(enrichment_data)
}
'
cat(improved_function)

cat("\n=== TESTING STEPS ===\n")
cat("1. Run the debug_enrichment_extraction.R script first\n")
cat("2. Check the data types and values in the ora_results table\n") 
cat("3. Apply the appropriate fix based on the findings\n")
cat("4. Test with the updated .extract_enrichment_data function\n")

cat("\n=== ADDITIONAL DIAGNOSTIC QUERIES ===\n")

cat("Check p_adjusted data types:\n")
cat("SELECT DISTINCT typeof(p_adjusted) FROM ora_results;\n\n")

cat("Check p_adjusted value ranges:\n")
cat("SELECT analysis_id, MIN(p_adjusted) as min_p, MAX(p_adjusted) as max_p FROM ora_results GROUP BY analysis_id;\n\n")

cat("Check for problematic values:\n")
cat("SELECT * FROM ora_results WHERE p_adjusted IS NULL OR p_adjusted = '' OR typeof(p_adjusted) != 'real';\n\n")

cat("Manual threshold check:\n")
cat("SELECT analysis_id, COUNT(*) as total, \n")
cat("       SUM(CASE WHEN CAST(p_adjusted AS REAL) <= 0.1 THEN 1 ELSE 0 END) as significant\n")
cat("FROM ora_results WHERE analysis_id IN (1, 2, 3) GROUP BY analysis_id;\n")

cat("\n=== END ANALYSIS ===\n")