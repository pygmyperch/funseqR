# GO Category Systems in funseqR

This document clarifies the two GO category formats used throughout the funseqR codebase.

## Two Category Systems

### 1. Database Storage Format (Single Letters)
- **Biological Process**: "P"
- **Molecular Function**: "F" 
- **Cellular Component**: "C"

**Used in:**
- Database `go_terms.go_category` column
- `process_annotations()` output `go_categories` field
- Internal data processing functions

**Source:** Extracted from GO term names as first character (line 1199 in `db_annotation.R`)

### 2. Analysis Interface Format (Full Names)
- **Biological Process**: "BP"
- **Molecular Function**: "MF"
- **Cellular Component**: "CC"

**Used in:**
- User-facing function parameters (e.g., `run_go_enrichment_analysis()`)
- Analysis result ontology labels
- External API compatibility

## Mapping Functions

### Consistent Mapping Pattern
```r
ontology_map <- c("BP" = "P", "MF" = "F", "CC" = "C")
category_code <- ontology_map[ontology]
```

### Functions with Proper Mapping
- `perform_go_enrichment()` (line 362 in `go_enrichment.R`) ✅
- `.perform_go_enrichment_from_data()` (line 436 in `enrichment_analysis.R`) ✅ **FIXED**

## Issue Resolved

**Problem:** `run_go_enrichment_analysis()` was looking for "BP"/"MF"/"CC" but database contained "P"/"F"/"C"

**Solution:** Added ontology mapping in `.perform_go_enrichment_from_data()` to convert analysis interface format to database storage format

**Commit:** 831f31e - Fix GO category mapping in enrichment analysis