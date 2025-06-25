# Coordinate Format Standardization Implementation

## Summary

Successfully implemented comprehensive coordinate format standardization across the funseqR package to support both `chrom` and `chromosome` column names while defaulting to the standard BED format (`chrom`).

## Changes Made

### 1. **vcf2bed_db() Function** (R/db_vcf.R:291-296)
- **Status**: ✅ Already correct
- **Output format**: `chrom, start, end` (standard BED format)
- **Implementation**: Creates BED data with `chrom` column name and 0-based start coordinates

### 2. **Enrichment Analysis Functions** (R/enrichment_analysis.R)

#### **run_go_enrichment_analysis()** (lines 38-163)
- **Status**: ✅ Updated with flexible coordinate detection
- **Implementation**: 
  - Detects both `chrom` and `chromosome` column names
  - Supports BED format: `chrom`/`chromosome`, `start`, `end`
  - Supports VCF format: `chrom`/`chromosome`, `position`
  - Provides helpful error messages showing actual column names found
  - Updated documentation to reflect flexible input formats

#### **run_kegg_enrichment_analysis()** (lines 197-297)
- **Status**: ✅ Updated with flexible coordinate detection
- **Implementation**: 
  - Same flexible coordinate detection as GO function
  - Supports both BED and VCF input formats
  - Updated documentation to reflect flexible input formats

### 3. **Documentation Updates**
- **Status**: ✅ Updated parameter descriptions
- **Changes**: Updated `@param candidate_loci` descriptions to specify support for both `chrom` and `chromosome` column names

## Implementation Details

### Coordinate Detection Logic
```r
# Detect available column names (support both standard formats)
has_chrom <- "chrom" %in% colnames(candidate_loci)
has_chromosome <- "chromosome" %in% colnames(candidate_loci)
has_start <- "start" %in% colnames(candidate_loci)
has_end <- "end" %in% colnames(candidate_loci)
has_position <- "position" %in% colnames(candidate_loci)

# Handle BED file format (chrom/chromosome, start, end)
if ((has_chrom || has_chromosome) && has_start && has_end) {
  chrom_col <- if (has_chrom) "chrom" else "chromosome"
  if (verbose) message("  - Detected BED format input (", chrom_col, ", start, end)")
  # ... processing logic
}
```

### Supported Input Formats

1. **BED format**: `chrom, start, end` (preferred) or `chromosome, start, end`
2. **VCF format**: `chrom, position` (preferred) or `chromosome, position`
3. **Locus ID vector**: Direct character vector of locus identifiers

### Error Handling
- Provides informative error messages showing actual column names found
- Suggests expected format options when input format is not recognized
- Maintains backward compatibility with existing `chromosome` column names

## Workflow Integration

### Data Flow
1. **VCF Import** → `import_vcf_to_db()` → stores with `chromosome, position`
2. **BED Creation** → `vcf2bed_db()` → creates `chrom, start, end` (standard format)
3. **Enrichment Analysis** → flexible functions accept both formats seamlessly

### Backward Compatibility
- ✅ Existing code using `chromosome` column names continues to work
- ✅ New code benefits from standard BED format (`chrom`)
- ✅ User can provide data in either format without modification

## Testing Scenarios Covered

1. **Standard BED input**: `data.frame(chrom = c("LG1"), start = c(1000), end = c(1001))`
2. **Legacy VCF input**: `data.frame(chromosome = c("LG1"), position = c(1001))`
3. **Mixed user input**: Functions detect and handle appropriately
4. **Error cases**: Clear error messages for unsupported formats

## Benefits Achieved

1. **Standardization**: Default to BED format standard (`chrom`) across the package
2. **Flexibility**: Support both formats for user convenience
3. **Clarity**: Verbose messages indicate which format was detected
4. **Future-proofing**: Easy to extend for additional coordinate formats
5. **Backward compatibility**: No breaking changes to existing workflows

## Files Modified

- `R/enrichment_analysis.R`: Updated both enrichment functions with flexible coordinate handling
- `R/db_vcf.R`: Already correctly using standard BED format
- Documentation: Updated parameter descriptions

## Status: ✅ Complete

The coordinate format standardization is fully implemented and tested. The package now supports both `chrom` and `chromosome` column names with intelligent detection, while defaulting to the standard BED format (`chrom, start, end`) for all new data creation.