#!/usr/bin/env Rscript

# Debug script to investigate enrichment data extraction mismatch
# Issue: ORA finds 2 significant terms but enrichment retrieval finds 0

library(DBI)
library(dplyr)

# Connect to database (adjust path as needed)
con <- dbConnect(RSQLite::SQLite(), "your_database.db")  # Replace with actual path

cat("=== Debugging Enrichment Data Extraction ===\n")

# Step 1: Check what's in the ORA tables
cat("\n1. Checking ORA analyses table:\n")
ora_analyses <- dbGetQuery(con, "SELECT * FROM ora_analyses ORDER BY analysis_id")
print(ora_analyses)

cat("\n2. Checking ORA results table:\n")
ora_results <- dbGetQuery(con, "SELECT * FROM ora_results ORDER BY analysis_id, p_adjusted")
print(ora_results)

# Step 3: Test the .extract_enrichment_data query with the exact parameters
cat("\n3. Testing enrichment data extraction with analysis_ids = c(1, 2, 3) and threshold = 0.1\n")

analysis_ids <- c(1, 2, 3)
significance_threshold <- 0.1

# Replicate the exact query from .extract_enrichment_data
analysis_ids_str <- paste(analysis_ids, collapse = ", ")
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
    AND res.p_adjusted <= ?
  ORDER BY res.p_adjusted
")

cat("Query being executed:\n")
cat(enrichment_query)
cat("\nParameters: significance_threshold =", significance_threshold, "\n")

enrichment_data <- dbGetQuery(con, enrichment_query, list(significance_threshold))

cat("\nResults from enrichment query:\n")
print(enrichment_data)
cat("Number of rows returned:", nrow(enrichment_data), "\n")

# Step 4: Let's check what we get without the p_adjusted filter
cat("\n4. Checking ALL results for these analysis IDs (no p_adjusted filter):\n")
all_results_query <- paste0("
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
  ORDER BY res.p_adjusted
")

all_results <- dbGetQuery(con, all_results_query)
print(all_results)

# Step 5: Check data types of p_adjusted column
cat("\n5. Checking data types and specific p_adjusted values:\n")
type_check_query <- paste0("
  SELECT 
    analysis_id,
    term_id,
    p_adjusted,
    typeof(p_adjusted) as p_adj_type,
    CASE 
      WHEN p_adjusted <= 0.1 THEN 'YES'
      ELSE 'NO'
    END as passes_threshold
  FROM ora_results 
  WHERE analysis_id IN (", analysis_ids_str, ")
  ORDER BY analysis_id, p_adjusted
")

type_check <- dbGetQuery(con, type_check_query)
print(type_check)

# Step 6: Let's manually check the comparison
cat("\n6. Manual comparison test:\n")
if (nrow(all_results) > 0) {
  for (i in 1:nrow(all_results)) {
    p_adj <- all_results$fdr[i]
    passes <- p_adj <= significance_threshold
    cat("Row", i, ": p_adjusted =", p_adj, ", <= 0.1 ?", passes, "\n")
  }
}

# Step 7: Check if there are any non-numeric values
cat("\n7. Checking for non-numeric p_adjusted values:\n")
non_numeric_check <- dbGetQuery(con, paste0("
  SELECT analysis_id, term_id, p_adjusted 
  FROM ora_results 
  WHERE analysis_id IN (", analysis_ids_str, ")
    AND (p_adjusted IS NULL OR p_adjusted = '' OR typeof(p_adjusted) != 'real')
"))
cat("Non-numeric p_adjusted values:\n")
print(non_numeric_check)

# Step 8: Try a different comparison approach
cat("\n8. Testing alternative comparison methods:\n")

# Test with string comparison
string_test_query <- paste0("
  SELECT analysis_id, term_id, p_adjusted
  FROM ora_results 
  WHERE analysis_id IN (", analysis_ids_str, ")
    AND CAST(p_adjusted AS REAL) <= 0.1
  ORDER BY p_adjusted
")

string_test <- dbGetQuery(con, string_test_query)
cat("Results with CAST(p_adjusted AS REAL) <= 0.1:\n")
print(string_test)

# Test with exact comparison
exact_test_query <- paste0("
  SELECT analysis_id, term_id, p_adjusted, (p_adjusted - 0.1) as difference
  FROM ora_results 
  WHERE analysis_id IN (", analysis_ids_str, ")
  ORDER BY p_adjusted
")

exact_test <- dbGetQuery(con, exact_test_query)
cat("Exact comparison test (difference from 0.1):\n")
print(exact_test)

# Clean up
dbDisconnect(con)

cat("\n=== Debug Complete ===\n")
cat("Summary of findings:\n")
cat("- ORA analyses found:", nrow(ora_analyses), "analyses\n")
cat("- ORA results found:", nrow(ora_results), "total results\n")
cat("- Enrichment query returned:", nrow(enrichment_data), "results\n")
cat("- Check the data types and comparison logic above to identify the issue\n")