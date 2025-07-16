# Critical ORA ClusterProfiler Fixes

## Root Causes Identified

The different ORA enrichment results between the cp_enrich and simple databases are caused by **TWO critical bugs** in the clusterProfiler implementation.

### Problem 1: Incorrect Universe Parameter

**Current (INCORRECT) implementation** in `.perform_clusterprofiler_enrichment()`:
```r
enrichment_result <- clusterProfiler::enricher(
  gene = go_data$foreground$genes,
  universe = go_data$background$genes,  # ❌ WRONG: excludes foreground genes
  TERM2GENE = clusterprofiler_data$term2gene,
  TERM2NAME = clusterprofiler_data$term2name,
  ...
)
```

### Problem 2: Incomplete TERM2GENE Mapping

**Current (INCORRECT) TERM2GENE building** in `.convert_go_data_to_clusterprofiler()`:
```r
for (gene in go_data$background$genes) {  # ❌ WRONG: only background genes
  if (gene %in% names(go_data$background$gene2go)) {
    gene_terms <- go_data$background$gene2go[[gene]]
    // ... build term2gene mapping
  }
}
```

**BOTH issues must be fixed** for clusterProfiler to work correctly:

### Why This Matters

In Over-Representation Analysis (ORA):
- **Foreground**: Genes of interest (candidate loci with annotations)  
- **Background**: Reference genes (non-candidate loci with annotations)
- **Universe**: ALL genes considered in the analysis (foreground + background)

The hypergeometric test calculates: "What's the probability of seeing X or more foreground genes in a term, given the term contains Y genes total in the universe?"

**With incorrect universe (background only):**
- Universe size is artificially reduced
- Statistical calculations are wrong
- P-values and enrichment significance change
- Different terms appear significant

### Comparison with Legacy Method

**Legacy method (CORRECT)** in `.perform_legacy_enrichment()`:
```r
# Uses separate counts for proper hypergeometric test
total_fg <- length(go_data$foreground$genes)
total_bg <- length(go_data$background$genes)

# Hypergeometric test with correct population sizes
p_value <- phyper(fg_with_term - 1, bg_with_term, total_bg - bg_with_term, total_fg, lower.tail = FALSE)
```

This correctly accounts for the full population.

### The Fixes

**Fix 1: Correct universe parameter:**
```r
enrichment_result <- clusterProfiler::enricher(
  gene = go_data$foreground$genes,
  universe = c(go_data$foreground$genes, go_data$background$genes),  # ✅ CORRECT: includes all genes
  TERM2GENE = clusterprofiler_data$term2gene,
  TERM2NAME = clusterprofiler_data$term2name,
  ...
)
```

**Fix 2: Complete TERM2GENE mapping:**
```r
# Include ALL genes (foreground + background) in TERM2GENE mapping
all_genes <- c(go_data$foreground$genes, go_data$background$genes)
all_gene2go <- c(go_data$foreground$gene2go, go_data$background$gene2go)

for (gene in all_genes) {  # ✅ CORRECT: all genes
  if (gene %in% names(all_gene2go)) {
    gene_terms <- all_gene2go[[gene]]
    // ... build complete term2gene mapping
  }
}
```

### Impact

This explains why:
1. **Different loci appear enriched** between versions
2. **Different p-values and FDR values** for the same terms  
3. **Only 3 out of 6 loci were common** between result sets

The cp_enrich database likely used the legacy method (correct statistics), while the simple database used the incorrect clusterProfiler universe.

### Next Steps

1. ✅ **Fix the universe parameter** in `.perform_clusterprofiler_enrichment()`
2. ✅ **Fix the TERM2GENE mapping** in `.convert_go_data_to_clusterprofiler()`
3. ✅ **Apply same fixes** to `.perform_clusterprofiler_kegg_enrichment()` and `.convert_kegg_data_to_clusterprofiler()`  
4. **Test the fixes** by re-running ORA and comparing results
5. **Verify** that fixed results match the legacy method results

These are critical bugs that completely invalidated the statistical validity of all clusterProfiler-based enrichment analyses. The TERM2GENE mapping bug was especially severe as it prevented clusterProfiler from seeing foreground genes at all.