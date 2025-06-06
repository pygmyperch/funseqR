# GO Enrichment Analysis in funseqR

This document describes the new GO enrichment functionality added to funseqR, which enables comparative analysis between candidate loci and background datasets to identify overrepresented biological processes, molecular functions, and cellular components.

## 🎯 Overview

The GO enrichment workflow is designed for scenarios where you have:
- **Background dataset**: Full set of variants (e.g., all SNPs from genome-wide analysis)
- **Foreground dataset**: Candidate variants of interest (e.g., adaptive loci, FST outliers, GWAS hits)

The analysis identifies which Gene Ontology terms are significantly overrepresented in your candidate loci compared to the background, helping to understand the biological processes potentially under selection or associated with your phenotype of interest.

## 🚀 Quick Start

### Simple Workflow

```r
library(funseqR)

# Connect to your annotated database
con <- connect_funseq_db("your_analysis.db")

# Run complete GO enrichment analysis
results <- run_go_enrichment_workflow(
  con = con,
  project_id = 1,
  candidate_vcf_file = "candidate_loci.vcf",
  background_file_id = 1,  # File ID of your full dataset
  ontologies = c("BP", "MF"),  # Biological Process + Molecular Function
  store_results = TRUE,
  create_plots = TRUE
)

# View results
print(results$summary)
print(results$plots$BP_bubble)  # Bubble plot
print(results$plots$BP_table)   # Summary table

close_funseq_db(con)
```

### Example: Snapper Adaptive Loci

```r
# Real example using snapper seascape genomics data
con <- connect_funseq_db("snapper_analysis.db")

results <- run_go_enrichment_workflow(
  con = con,
  project_id = 1,
  candidate_vcf_file = "SA448_855.vcf",     # 855 candidate adaptive SNPs
  background_file_id = 1,                    # SA448_14699.vcf (full 14,699 SNPs)
  ontologies = c("BP", "MF"),
  min_genes = 5,
  verbose = TRUE
)

# Results show enrichment for processes like:
# - response to temperature stimulus
# - ion transport
# - metabolic processes
# - enzyme activity
```

## 📊 Visualization Options

### 1. Bubble Plots
Show GO terms ranked by fold enrichment, with bubble size = gene count and color = significance:

```r
bubble_plot <- create_go_bubble_plot(
  enrichment_results = bp_results,
  max_terms = 20,
  min_fold_enrichment = 1.5
)
```

### 2. Treemaps
Hierarchical view where area represents gene count:

```r
# Requires treemapify package
treemap_plot <- create_go_treemap(bp_results)
```

### 3. Multi-Ontology Comparisons
Compare across Biological Process, Molecular Function, and Cellular Component:

```r
comparison_plot <- create_go_comparison_plot(bp_results, mf_results, cc_results)
```

### 4. Interactive Plots
Hover-enabled plots with plotly:

```r
# Requires plotly package
interactive_plot <- create_interactive_go_plot(bp_results, max_terms = 25)
```

## 🔧 Core Functions

### High-Level Workflow
- `run_go_enrichment_workflow()`: Complete end-to-end analysis
- `generate_go_enrichment_report_section()`: Generate R Markdown for reports

### Step-by-Step Functions
- `import_candidate_loci()`: Import and link candidate VCF to annotations
- `extract_go_terms_for_enrichment()`: Get GO terms for foreground/background
- `perform_go_enrichment()`: Statistical testing with hypergeometric tests
- `store_go_enrichment_results()`: Save results to database
- `get_go_enrichment_results()`: Retrieve stored results

### Visualization Functions  
- `create_go_bubble_plot()`: Standard bubble plot visualization
- `create_go_treemap()`: Treemap visualization
- `create_go_summary_table()`: Publication-ready summary tables
- `create_go_comparison_plot()`: Multi-ontology comparison
- `create_interactive_go_plot()`: Interactive plotly visualization

## 📊 Statistical Methods

### Hypergeometric Enrichment Test
For each GO term, tests whether it's overrepresented in foreground vs background:

- **Null hypothesis**: GO term frequency is the same in foreground and background
- **Test statistic**: Number of foreground genes with the GO term
- **P-value**: Hypergeometric probability P(X ≥ observed | population, successes, sample_size)
- **Multiple testing**: FDR correction (Benjamini-Hochberg)

### Parameters
- `min_genes`: Minimum genes required for testing (default: 5)
- `max_genes`: Maximum genes to exclude very broad terms (default: 500)
- `ontology`: "BP" (Biological Process), "MF" (Molecular Function), "CC" (Cellular Component)

### Interpretation
- **Fold Enrichment**: Observed/Expected ratio
  - 1.5-2x: Modest enrichment
  - 2-5x: Strong enrichment  
  - >5x: Very strong enrichment
- **Significance**: FDR < 0.05 (significant), FDR < 0.01 (highly significant)

## 🗃️ Database Schema

New tables for storing enrichment results:

```sql
-- Analysis metadata
CREATE TABLE go_enrichment_analyses (
  enrichment_id INTEGER PRIMARY KEY,
  project_id INTEGER NOT NULL,
  foreground_file_id INTEGER NOT NULL,
  background_file_id INTEGER NOT NULL,
  ontology TEXT NOT NULL,
  analysis_date TEXT NOT NULL,
  total_foreground_genes INTEGER,
  total_background_genes INTEGER,
  analysis_parameters TEXT
);

-- Individual term results
CREATE TABLE go_enrichment_results (
  result_id INTEGER PRIMARY KEY,
  enrichment_id INTEGER NOT NULL,
  go_id TEXT NOT NULL,
  go_term TEXT NOT NULL,
  go_category TEXT NOT NULL,
  foreground_count INTEGER NOT NULL,
  background_count INTEGER NOT NULL,
  total_foreground INTEGER NOT NULL,
  total_background INTEGER NOT NULL,
  expected_count REAL,
  fold_enrichment REAL,
  p_value REAL NOT NULL,
  p_adjusted REAL,
  significance_level TEXT
);
```

## 📚 Integration with Analysis Reports

Add GO enrichment sections to your existing funseqR reports:

```r
# Automatically detects candidate/background files and generates plots
go_section <- generate_go_enrichment_report_section(
  con = con,
  project_id = 1,
  candidate_file_pattern = "candidate|adaptive|outlier",
  max_terms_plot = 15,
  include_treemap = TRUE
)

# Returns R Markdown code for inclusion in reports
```

## 🔍 Use Cases

### 1. Adaptive Evolution Studies
- **Foreground**: FST outliers, environmental association SNPs
- **Background**: Genome-wide SNPs
- **Goal**: Identify adaptive biological processes

### 2. GWAS Follow-up
- **Foreground**: Significant GWAS hits
- **Background**: All tested SNPs
- **Goal**: Understand biological basis of trait associations

### 3. Population Genomics
- **Foreground**: Differentiated loci between populations
- **Background**: All variable sites
- **Goal**: Characterize divergent selection pressures

### 4. Conservation Genomics
- **Foreground**: Loci associated with fitness/adaptation
- **Background**: Neutral variation
- **Goal**: Prioritize genes for conservation

## 📦 Dependencies

### Required
- `DBI`, `RSQLite`: Database operations
- `dplyr`: Data manipulation
- `ggplot2`: Visualization
- `stringr`: String operations

### Optional (for enhanced features)
- `treemapify`: Treemap visualizations
- `plotly`: Interactive plots  
- `knitr`, `rmarkdown`: Report generation

## 🛠️ Installation & Setup

The GO enrichment functionality is included in the main funseqR package. To use it:

1. **Complete a standard funseqR analysis** with VCF import, BLAST, and annotation
2. **Prepare your candidate dataset** (subset VCF of interesting loci)
3. **Run GO enrichment** using the workflow functions

```r
# Install development version with GO enrichment
devtools::install_github("user/funseqR", ref = "go-enrichment-feature")
```

## 📋 Troubleshooting

### Common Issues

**No significant results**
- Check that background and foreground datasets are appropriate
- Ensure sufficient annotation coverage
- Consider relaxing significance thresholds

**Very broad GO terms dominate**
- Increase `max_genes` parameter to exclude broad terms
- Focus on more specific ontology levels

**Low gene overlap**
- Verify that candidate and background use same reference genome
- Check annotation quality and coverage
- Consider expanding flanking regions for BLAST

**High fold enrichments (>10x)**
- May indicate technical artifacts
- Check for batch effects or annotation biases
- Validate with independent datasets

### Getting Help

1. Check the vignette: `vignette("go_enrichment_analysis", package = "funseqR")`
2. View function documentation: `?run_go_enrichment_workflow`
3. Run the example script: `examples/go_enrichment_example.R`

## 🎓 Citation

If you use the GO enrichment functionality in your research, please cite funseqR and mention the specific features used:

```
Used funseqR GO enrichment analysis to identify overrepresented biological 
processes in candidate adaptive loci using hypergeometric testing with FDR 
correction (Benjamini-Hochberg).
```

## 🔮 Future Enhancements

Planned features for future releases:
- KEGG pathway enrichment analysis  
- Protein-protein interaction networks
- Gene set enrichment analysis (GSEA)
- Integration with external enrichment databases
- Automated report generation with publication-ready figures

---

**Note**: This is a development feature currently on the `go-enrichment-feature` branch. It will be merged to main after thorough testing and validation.