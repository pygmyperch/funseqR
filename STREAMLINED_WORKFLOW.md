# FunseqR Streamlined Workflow

This document shows the updated streamlined workflow with the cleaned-up function names.

## Updated Workflow Example

```r
library(funseqR)

# Define input paths
db_path <- "funseq_project.db"
vcf_file <- "SA448_14699.vcf"
ref_genome_file <- "Chrysophrys_auratus.v.1.0.all.assembly.units.fna"
blast_db_path <- "/Volumes/SSD2TB/blastDBs/teleost"
blast_db_name <- "teleostei_db"

# 1. Create database 
con <- create_funseq_db(db_path, force = TRUE)

# 2. Import data
vcf_import <- import_vcf(con, vcf_file)
ref_import <- import_reference(con, ref_genome_file, "Chrysophrys_auratus", "v1.0")

# 3. Import locus statistics and define candidate loci
define_locus_statistics(con,
                        statistics = my_pvalues_df,
                        candidate_threshold = 0.01)

# 4. Clean unmapped scaffolds
main_chroms <- c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", 
                 "LG9", "LG10", "LG11", "LG12", "LG13", "LG14", "LG15", "LG16", 
                 "LG17", "LG18", "LG19", "LG20", "LG21", "LG22", "LG23", "LG24")
define_chromosomes(con, main_chroms)

# 5. Extract flanking sequences
flanking_import <- import_flanking_seqs(con,
                                        vcf_import$file_id,
                                        ref_import$genome_id,
                                        flank_size = 300,
                                        translate_flanks = TRUE,
                                        threads = 6,
                                        batch_size = 1000,
                                        verbose = TRUE)

# 6. Perform BLAST
blast_results <- perform_blast(con,
                               vcf_import$file_id,
                               db_path = blast_db_path,
                               db_name = blast_db_name,
                               blast_type = "diamond_blastx",
                               seq_type = "orf_nuc",
                               e_value = 1e-10,
                               max_hits = 5,
                               threads = 4)

# 7. Annotate BLAST results
annotation_results <- annotate_blast_results(con,
                                             blast_results$blast_param_id,
                                             evidence_keep = c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "TAS", "IC", "IEA", "ISS"),
                                             max_hits = 1,
                                             verbose = TRUE)

# 8. Generate base annotation data
results <- compile_funseq_results(
  con, 
  stage = "annotations",
  include = c("GO", "KEGG", "Pfam", "InterPro", "eggNOG"),
  export_csv = "step1_annotations.csv"
)

# 9. Run ORA enrichment
ORA_results <- run_ORA(con,
                       annotation_type = "all",
                       significance_threshold = 0.1,
                       blast_param_id = 1)

# 10. Add enrichment results
results <- compile_funseq_results(
  con,
  stage = "enrichment", 
  data = results,
  analysis_ids = c(1, 2, 3),
  significance_threshold = 0.1,
  export_csv = "step2_with_enrichment.csv",
  verbose = TRUE
)

# 11. Create Manhattan plot (using stored statistics automatically)
manhattan_plot <- create_functional_manhattan_plot(
  con,
  # y_values automatically retrieved from stored statistics
  y_label = "q-value",
  vcf_file_id = 1,
  enrichment_data = results,
  label_type = "gene_name",
  label_cex = 1.5,
  numeric_x_labels = TRUE,
  enriched_point_size = 2,
  enriched_point_shape = 17,  # triangle
  enriched_point_color = "red",
  use_label_lines = TRUE
)

print(manhattan_plot)
ggsave("RDA_manhattan_streamlined.pdf", width = 8, height = 6)

# 12. Export GO terms for ReviGO (using stored candidates automatically)
export_revigo_file(con, 
                   source_type = "candidate_loci",
                   include_pvalues = FALSE,
                   output_file = "candidate_revigo.txt")
```

## Key Improvements

### 🎯 **Simplified Function Names**
- `import_vcf_to_db()` → `import_vcf()`
- `import_reference_to_db()` → `import_reference()`
- `import_flanking_seqs_to_db()` → `import_flanking_seqs()`
- `perform_blast_db()` → `perform_blast()`
- `vcf2bed_db()` → `vcf_to_bed()`
- `create_funseq_schema()` → `create_schema()`

### 🗂️ **Cleaner File Structure**
- Removed 12 obsolete files (~50% reduction)
- Renamed 7 files (removed `db_` prefix)
- 14 focused files remaining

### ⚡ **Streamlined Workflow**
- **Stored Statistics**: `y_values` automatically retrieved from `define_locus_statistics()`
- **Stored Candidates**: Automatically used in `export_revigo_file()` 
- **One Database = One Project**: Clean separation of analyses
- **Intelligent Defaults**: Functions automatically use stored data when available

### 📊 **Package Stats**
- **Files**: 14 (down from 26)
- **Exported Functions**: 56 (focused on essential workflow)
- **Core Workflow Functions**: All 12 key functions available ✅

The package is now streamlined, focused, and maintains all functionality while being much easier to use and understand!