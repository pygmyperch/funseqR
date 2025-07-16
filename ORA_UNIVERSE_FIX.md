# Critical ORA Universe Fix

## Root Cause Identified

The different ORA enrichment results between the cp_enrich and simple databases are caused by an **incorrect universe parameter** in the clusterProfiler implementation.

### The Problem

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

**The universe should include ALL genes** (both foreground and background) for proper statistical testing.

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

### The Fix

**Correct clusterProfiler implementation:**
```r
enrichment_result <- clusterProfiler::enricher(
  gene = go_data$foreground$genes,
  universe = c(go_data$foreground$genes, go_data$background$genes),  # ✅ CORRECT: includes all genes
  TERM2GENE = clusterprofiler_data$term2gene,
  TERM2NAME = clusterprofiler_data$term2name,
  ...
)
```

### Impact

This explains why:
1. **Different loci appear enriched** between versions
2. **Different p-values and FDR values** for the same terms  
3. **Only 3 out of 6 loci were common** between result sets

The cp_enrich database likely used the legacy method (correct statistics), while the simple database used the incorrect clusterProfiler universe.

### Next Steps

1. **Fix the universe parameter** in `.perform_clusterprofiler_enrichment()`
2. **Apply same fix** to `.perform_clusterprofiler_kegg_enrichment()`  
3. **Test the fix** by re-running ORA and comparing results
4. **Verify** that fixed results match the legacy method results

This is a critical bug that affects the statistical validity of all clusterProfiler-based enrichment analyses.