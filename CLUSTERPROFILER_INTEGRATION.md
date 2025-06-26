# clusterProfiler Integration for funseqR

This document describes the new clusterProfiler integration for GO enrichment analysis in funseqR.

## Overview

The clusterProfiler integration provides a robust, well-tested alternative to the custom enrichment functions in funseqR. It uses the same statistical methodology (hypergeometric test) as the original funseqR functions but leverages clusterProfiler's mature implementation and enhanced visualization capabilities.

## Key Benefits

1. **Proven Statistical Foundation**: Uses the same hypergeometric test as original funseqR
2. **Enhanced Reliability**: Battle-tested implementation with proper edge case handling
3. **Better Visualizations**: Built-in plots (dotplot, barplot, etc.) with publication-quality output
4. **Consistent Results**: Should reproduce the original results (29 BP terms tested, 2 significant at FDR < 0.1)
5. **Future-Proof**: Actively maintained by the Bioconductor community

## Installation Requirements

```r
# Install required Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "DOSE"))
```

## Quick Start

After completing your standard funseqR workflow (VCF import, BLAST, annotation):

```r
# Run clusterProfiler enrichment analysis
enrichment_results <- run_clusterprofiler_enrichment(
  con = con,
  candidate_file_id = vcf_cand_import$file_id,
  background_file_id = vcf_import$file_id,
  blast_param_id = blast_results$blast_param_id,
  ontologies = c("BP", "MF", "CC"),
  pvalue_cutoff = 0.1  # To match previous 0.1 threshold
)

# View results
print(enrichment_results)

# Create visualizations
if (!is.null(enrichment_results$BP)) {
  dotplot(enrichment_results$BP)
  barplot(enrichment_results$BP)
}
```

## Function Reference

### `run_clusterprofiler_enrichment()`

Main function for running GO enrichment analysis using clusterProfiler.

**Parameters:**
- `con`: Database connection object
- `candidate_file_id`: File ID of candidate dataset (e.g., adaptive loci)
- `background_file_id`: File ID of background dataset (e.g., all sequenced loci)
- `blast_param_id`: BLAST parameter ID for annotations
- `ontologies`: GO ontologies to test ("BP", "MF", "CC")
- `pvalue_cutoff`: P-value cutoff for significance (default: 0.05)
- `padjust_method`: P-value adjustment method (default: "BH")
- `min_gs_size`: Minimum gene set size (default: 5)
- `max_gs_size`: Maximum gene set size (default: 500)

**Returns:**
List containing clusterProfiler enrichResult objects for each ontology plus summary information.

### `compare_enrichment_methods()`

Compare clusterProfiler results with original funseqR workflow for validation.

### `convert_clusterprofiler_to_funseqr()`

Convert clusterProfiler results to original funseqR format for backward compatibility.

### `generate_enrichment_report()`

Generate HTML report with visualizations and tables.

## Expected Results

Based on the original funseqR workflow, the clusterProfiler integration should produce:

- **BP terms tested**: ~29
- **Significant BP terms (FDR < 0.1)**: 2
- **Background genes**: 527
- **Foreground genes**: 39

## Workflow Integration

The clusterProfiler integration fits seamlessly into the existing funseqR workflow:

```r
# 1. Standard funseqR setup (unchanged)
con <- create_funseq_db(db_path, force = TRUE)
vcf_import <- import_vcf_to_db(con, vcf_file)
vcf_cand_import <- import_vcf_to_db(con, vcf_file_cand)
# ... (rest of standard workflow)

# 2. BLAST and annotation (unchanged)
blast_results <- perform_blast_db(con, ...)
annotation_results <- annotate_blast_results(con, ...)

# 3. NEW: clusterProfiler enrichment
enrichment_results <- run_clusterprofiler_enrichment(
  con = con,
  candidate_file_id = vcf_cand_import$file_id,
  background_file_id = vcf_import$file_id,
  blast_param_id = blast_results$blast_param_id
)

# 4. Analysis and visualization
print(enrichment_results)
dotplot(enrichment_results$BP)
```

## Visualization Options

clusterProfiler provides many built-in visualization functions:

```r
# Dot plot (shows gene ratio and p-values)
dotplot(enrichment_results$BP, showCategory = 15)

# Bar plot (shows gene counts)
barplot(enrichment_results$BP, showCategory = 10)

# Network plot (requires additional setup)
# emapplot(enrichment_results$BP)

# Gene-concept network
# cnetplot(enrichment_results$BP)
```

## Comparison with Original Functions

| Feature | Original funseqR | clusterProfiler |
|---------|------------------|-----------------|
| Statistical test | Hypergeometric | Hypergeometric (same) |
| Gene identification | UniProt accession | UniProt accession (same) |
| Multiple testing | FDR correction | FDR correction (same) |
| Visualization | Custom plots | Rich built-in plots |
| Maintenance | Custom code | Community maintained |
| Edge cases | Manual handling | Robust handling |

## Files Added

- `R/clusterprofiler_enrichment.R`: Main integration functions
- `R/enrichment_comparison.R`: Comparison and utility functions  
- `examples/clusterprofiler_workflow_example.R`: Complete workflow example
- `test_clusterprofiler_integration.R`: Test script

## Troubleshooting

### Package Installation Issues
```r
# If installation fails, try:
BiocManager::install("clusterProfiler", force = TRUE)
BiocManager::install("DOSE", force = TRUE)
```

### No Significant Results
- Check p-value cutoff (try 0.1 instead of 0.05)
- Verify gene lists are properly extracted
- Check that GO annotations exist in database

### Different Results from Original
- Compare gene counts using `compare_enrichment_methods()`
- Verify same blast_param_id is used
- Check ontology filtering

## Migration from Original Functions

To migrate from the original `run_go_enrichment_workflow()`:

1. Replace the workflow call with `run_clusterprofiler_enrichment()`
2. Update visualization code to use clusterProfiler functions
3. Adjust p-value thresholds if needed
4. Use comparison functions to validate results

## Support

For questions about the clusterProfiler integration, check:
1. The example workflow in `examples/clusterprofiler_workflow_example.R`
2. The test script `test_clusterprofiler_integration.R`
3. clusterProfiler documentation: https://bioconductor.org/packages/clusterProfiler/
4. funseqR documentation and vignettes