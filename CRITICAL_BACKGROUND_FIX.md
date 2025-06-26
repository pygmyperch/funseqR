# Critical Background Set Fix for Enrichment Analysis

## Problem Identified

The new enrichment functions were using **incorrect background sets**, causing statistical tests to fail detecting real enrichment.

### Root Cause
- **New functions**: Used all 694 annotations as background (including those without GO/KEGG terms)
- **Old functions**: Used only annotations with GO/KEGG terms as background (527 for GO)
- **Impact**: Statistical tests were heavily diluted, making enrichment undetectable

## Evidence of the Problem

| Function | Background Count | Terms Tested | Significant Results |
|----------|------------------|--------------|-------------------|
| **Old (working)** | 527 (GO only) | 29 BP terms | **2 significant** |
| **New (broken)** | 694 (all annotations) | 13 BP terms | **0 significant** |

## Fixes Implemented

### 1. **Fixed Background Set Definition**

**Before (incorrect)**:
```r
background = list(
  loci = annotations$locus_id,  # All 694 annotations
  genes = annotations$locus_id
)
```

**After (correct)**:
```r
# Filter to only annotations with GO terms first
go_annotations <- annotations[!is.na(annotations$go_terms) & annotations$go_terms != "", ]

background = list(
  loci = go_annotations$locus_id,  # Only 527 with GO terms
  genes = go_annotations$locus_id
)
```

### 2. **Applied Same Fix to KEGG Function**

Applied identical logic to KEGG enrichment to maintain consistency.

### 3. **Added Debugging Output**

```r
if (verbose) {
  message("    - Total annotations: ", nrow(annotations))
  message("    - Annotations with GO terms: ", nrow(go_annotations))
}
```

### 4. **Cleaned Up Redundant Filtering**

Removed duplicate filtering in main functions since helper functions now handle it properly.

## Expected Results

With these fixes, the new enrichment functions should now:

1. **Match old function results**: Find the same 2 significant BP terms
2. **Use correct background**: 527 annotations with GO terms (not 694 total)
3. **Test more terms**: Should test ~29 BP terms (not just 13)
4. **Show proper statistics**: Realistic enrichment ratios and p-values

## Test Command

```r
# This should now find 2 significant BP terms at FDR < 0.1
go_results <- run_go_enrichment_analysis(
  annotations, 
  candidate_loci = candidates$bed_file,
  ontologies = c("BP", "MF", "CC"),
  significance_threshold = 0.1
)
```

## Files Modified

- **R/enrichment_analysis.R**:
  - Fixed `.prepare_go_data_from_annotations()` background definition
  - Fixed `.prepare_kegg_data_from_annotations()` background definition  
  - Added verbose debugging output
  - Removed redundant filtering in main functions

This fix addresses the fundamental statistical issue that was preventing enrichment detection.