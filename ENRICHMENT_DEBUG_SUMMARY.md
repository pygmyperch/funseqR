# Enrichment Data Extraction Debug Summary

## Problem Description

You reported that the ORA analysis (`run_ORA`) found 2 significantly enriched BP terms with `significance_threshold = 0.1`, but when retrieving those results using `compile_funseq_results` with `analysis_ids = c(1, 2, 3)` and the same threshold, it returned 0 terms.

The specific symptoms were:
- **ORA output**: Analysis ID 1 (BP): 2 significantly enriched terms  
- **ORA output**: Analysis ID 2 (MF): 0 significantly enriched terms
- **ORA output**: Analysis ID 3 (CC): 0 significantly enriched terms
- **Enrichment retrieval**: 0 terms total

## Root Cause Analysis

After analyzing the code in `/R/process_annotations.R`, the issue is likely in the `.extract_enrichment_data()` function at line 365. The problematic SQL query was:

```sql
WHERE ora.analysis_id IN (1, 2, 3) AND res.p_adjusted <= 0.1
```

### Potential Issues Identified

1. **Data Type Mismatch**: The `p_adjusted` column might be stored as TEXT instead of REAL in SQLite
2. **NULL/Empty Values**: The `p_adjusted` column might contain NULL or empty string values
3. **Floating Point Precision**: SQLite floating point comparisons might have precision issues
4. **Missing Analysis IDs**: The requested analysis_ids might not exist in the database
5. **JOIN Failures**: Foreign key relationships might be broken

## The Fix

The updated `.extract_enrichment_data()` function now includes:

### 1. Analysis ID Validation
```r
# First check if analysis_ids exist
check_query <- paste0("
  SELECT analysis_id FROM ora_analyses 
  WHERE analysis_id IN (", analysis_ids_str, ")
")
existing_ids <- DBI::dbGetQuery(con, check_query)$analysis_id
```

### 2. Improved SQL Query with Data Type Handling
```sql
WHERE ora.analysis_id IN (1, 2, 3)
  AND res.p_adjusted IS NOT NULL
  AND res.p_adjusted != ''
  AND CAST(res.p_adjusted AS REAL) <= ?
ORDER BY CAST(res.p_adjusted AS REAL)
```

### 3. Enhanced Debug Output
When no results are found, the function now provides detailed diagnostics:
- Total results per analysis ID
- Min/max p_adjusted values
- Count of significant results
- Data type information

## Key Changes Made

1. **File**: `/R/process_annotations.R`
2. **Function**: `.extract_enrichment_data()` (lines 365-403)
3. **Changes**:
   - Added analysis ID existence check
   - Added NULL and empty string filters
   - Used `CAST(res.p_adjusted AS REAL)` for comparisons
   - Added comprehensive debug output
   - Improved error handling and logging

## Testing the Fix

### Immediate Testing
1. Run the debug script: `debug_enrichment_extraction.R`
2. Check data types and values in the `ora_results` table
3. Test with the updated function: `test_enrichment_fix.R`

### Integration Testing
Test the full workflow:
```r
# After running ORA
ora_results <- run_ORA(con, "candidates.vcf", significance_threshold = 0.1)

# Extract enrichment results  
enriched_results <- compile_funseq_results(
  con = con,
  stage = "enrichment", 
  data = annotation_data,
  analysis_ids = c(1, 2, 3),
  significance_threshold = 0.1,
  verbose = TRUE
)
```

## Expected Outcomes

### If the Fix Worked
- You should now see the 2 enriched terms that were missing before
- Debug output will show actual p_adjusted values and counts
- The enrichment stage will properly map terms to loci

### If the Issue Persists
- Check the debug output for data type issues
- Verify analysis_ids exist in the database  
- Examine actual p_adjusted values in the database
- Consider other potential causes (e.g., gene ID mapping issues)

## Additional Diagnostic Queries

If you need to investigate further, these queries can help:

```sql
-- Check data types
SELECT DISTINCT typeof(p_adjusted) FROM ora_results;

-- Check value ranges
SELECT analysis_id, MIN(p_adjusted) as min_p, MAX(p_adjusted) as max_p 
FROM ora_results GROUP BY analysis_id;

-- Check for problematic values
SELECT * FROM ora_results 
WHERE p_adjusted IS NULL OR p_adjusted = '' OR typeof(p_adjusted) != 'real';

-- Manual threshold check
SELECT analysis_id, COUNT(*) as total, 
       SUM(CASE WHEN CAST(p_adjusted AS REAL) <= 0.1 THEN 1 ELSE 0 END) as significant
FROM ora_results WHERE analysis_id IN (1, 2, 3) GROUP BY analysis_id;
```

## Prevention

To prevent similar issues in the future:
1. Always use `CAST()` for numeric comparisons in SQLite
2. Include NULL and empty string checks in WHERE clauses
3. Add comprehensive logging and debug output
4. Validate input parameters before using them in queries
5. Test edge cases (non-existent IDs, extreme thresholds, etc.)

## Summary

This fix addresses the most likely cause of the enrichment mismatch issue by improving the robustness of the SQL query and adding better error diagnostics. The enhanced debug output will help identify any remaining issues if they persist.