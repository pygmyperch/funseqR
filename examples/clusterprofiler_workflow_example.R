# clusterProfiler Integration Example
# This script demonstrates the complete workflow using clusterProfiler for enrichment analysis

# Load required libraries
library(funseqR)
library(clusterProfiler)  # Make sure this is installed: BiocManager::install("clusterProfiler")
library(DOSE)             # Make sure this is installed: BiocManager::install("DOSE")

# === PART 1: Standard funseqR Workflow ===
# (This is your existing workflow - no changes needed)

# Define paths
db_path <- "funseq_project.db"
vcf_file <- "SA448_14699.vcf"          # Background/reference dataset
vcf_file_cand <- "SA448_855.vcf"       # Candidate/adaptive loci
ref_genome_file <- "Chrysophrys_auratus.v.1.0.all.assembly.units.fna"
blast_db_path <- "/Volumes/SSD2TB/blastDBs/teleost"
blast_db_name <- "teleostei_db"

# 1. Create database 
con <- create_funseq_db(db_path, force = TRUE)

# 2. Import data
vcf_import <- import_vcf_to_db(con, vcf_file)
vcf_cand_import <- import_vcf_to_db(con, vcf_file_cand)
ref_import <- import_reference_to_db(con, ref_genome_file, "Chrysophrys_auratus", "v1.0")

# 3. Define main chromosomes
main_chroms <- c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", 
                 "LG11", "LG12", "LG13", "LG14", "LG15", "LG16", "LG17", "LG18", "LG19", 
                 "LG20", "LG21", "LG22", "LG23", "LG24")
define_chromosomes(con, main_chroms)

# 4. Extract flanking sequences
flanking_import <- import_flanking_seqs_to_db(con,
                                              vcf_import$file_id,
                                              ref_import$genome_id,
                                              flank_size = 300,
                                              translate_flanks = TRUE,
                                              threads = 6,
                                              batch_size = 1000,
                                              verbose = TRUE)

# 5. Perform BLAST searches
blast_results <- perform_blast_db(con,
                                  vcf_import$file_id,
                                  db_path = blast_db_path,
                                  db_name = blast_db_name,
                                  blast_type = "diamond_blastx",
                                  seq_type = "orf_nuc",
                                  e_value = 1e-10,
                                  max_hits = 5,
                                  threads = 4)

# 6. Annotate BLAST results
annotation_results <- annotate_blast_results(con,
                                             blast_results$blast_param_id,
                                             evidence_keep = c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "TAS", "IC", "IEA", "ISS"),
                                             max_hits = 1,
                                             verbose = TRUE)

# === PART 2: NEW clusterProfiler Enrichment Analysis ===

# Run enrichment analysis using clusterProfiler
enrichment_results <- run_clusterprofiler_enrichment(
  con = con,
  candidate_file_id = vcf_cand_import$file_id,      # Candidate dataset
  background_file_id = vcf_import$file_id,          # Background dataset  
  blast_param_id = blast_results$blast_param_id,    # Annotation results
  ontologies = c("BP", "MF", "CC"),                 # All ontologies
  pvalue_cutoff = 0.1,                              # More lenient threshold to match old results
  padjust_method = "BH",                            # Benjamini-Hochberg correction
  min_gs_size = 5,                                  # Minimum 5 genes per term (matches old function)
  max_gs_size = 500,                                # Maximum 500 genes per term
  verbose = TRUE
)

# === PART 3: View and Analyze Results ===

# Print overall summary
print(enrichment_results)

# View BP (Biological Process) results
if (!is.null(enrichment_results$BP)) {
  cat("\\n=== Biological Process Results ===\\n")
  print(head(enrichment_results$BP@result))
  
  # Create visualizations
  if (nrow(enrichment_results$BP@result) > 0) {
    # Dot plot (shows top enriched terms)
    p1 <- dotplot(enrichment_results$BP, showCategory = 15)
    print(p1)
    
    # Bar plot (shows enrichment strength)
    p2 <- barplot(enrichment_results$BP, showCategory = 10)
    print(p2)
  }
}

# View MF (Molecular Function) results
if (!is.null(enrichment_results$MF)) {
  cat("\\n=== Molecular Function Results ===\\n")
  print(head(enrichment_results$MF@result))
}

# View CC (Cellular Component) results  
if (!is.null(enrichment_results$CC)) {
  cat("\\n=== Cellular Component Results ===\\n")
  print(head(enrichment_results$CC@result))
}

# === PART 4: Export Results ===

# Export significant results to CSV
if (!is.null(enrichment_results$BP) && nrow(enrichment_results$BP@result) > 0) {
  # Filter for significant results
  bp_sig <- enrichment_results$BP@result[enrichment_results$BP@result$p.adjust < 0.1, ]
  if (nrow(bp_sig) > 0) {
    write.csv(bp_sig, "BP_enrichment_results.csv", row.names = FALSE)
    cat("\\nSignificant BP results exported to BP_enrichment_results.csv\\n")
  }
}

# === PART 5: Compare with Old Results (Validation) ===

cat("\\n=== Validation Against Old Results ===\\n")
cat("Expected from old function:\\n")
cat("  - BP terms tested: 29\\n")
cat("  - Significant BP terms (FDR < 0.1): 2\\n")
cat("  - Background genes: 527\\n")
cat("  - Foreground genes: 39\\n\\n")

cat("New clusterProfiler results:\\n")
if (!is.null(enrichment_results$BP)) {
  bp_tested <- nrow(enrichment_results$BP@result)
  bp_significant <- sum(enrichment_results$BP@result$p.adjust < 0.1, na.rm = TRUE)
  cat("  - BP terms tested:", bp_tested, "\\n")
  cat("  - Significant BP terms (FDR < 0.1):", bp_significant, "\\n")
}
cat("  - Background genes:", enrichment_results$summary$background_genes, "\\n")
cat("  - Foreground genes:", enrichment_results$summary$candidate_genes, "\\n")

# Close database connection
dbDisconnect(con)

cat("\\n=== Analysis Complete ===\\n")
cat("Results available in 'enrichment_results' object\\n")
cat("Use dotplot(), barplot(), or other clusterProfiler visualization functions\\n")