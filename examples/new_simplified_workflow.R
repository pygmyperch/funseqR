# funseqR New Simplified Workflow Example
# =====================================

# This script demonstrates the new simplified workflow for functional annotation analysis
# The new approach dramatically reduces complexity while maintaining all analytical power

library(funseqR)

# ========================================
# STEP 1: BLAST and Annotation (unchanged)
# ========================================

# Run BLAST analysis (same as before)
blast_results <- run_blast(sequences, database)

# Annotate BLAST results (same as before)  
con <- annotate_blast_results(blast_results)

# ========================================
# STEP 2: Process Annotations (NEW - replaces many functions!)
# ========================================

# OLD WAY (complex - multiple functions needed):
# go_summary <- summarize_go_annotations(con, candidate_loci = candidates)
# kegg_summary <- summarize_kegg_pathways(con, candidate_loci = candidates, use_keggrest_enhancement = TRUE)
# profile <- generate_candidate_loci_profile(con, candidate_loci = candidates)
# # ... many more steps

# NEW WAY (simple - one function does it all):
annotations <- process_annotations(
  con, 
  include = c("GO", "KEGG"),                    # Default: both GO and KEGG
  export_csv = "my_functional_annotations.csv" # Optional: export for external tools
)

# The result is a comprehensive data frame with ALL functional annotations:
# - locus_id, chromosome, position (genomic coordinates)
# - gene_name, protein_name, uniprot_accession (protein info)
# - go_terms, go_names, go_categories (GO annotations)
# - kegg_pathways, kegg_pathway_names (KEGG with enhanced names using KEGGREST)
# - kegg_brite_categories (official KEGG BRITE classifications)
# - kegg_modules (functional module assignments)

# Preview the results
head(annotations)
summary(annotations)

# ========================================
# STEP 3: Enrichment Analysis (simplified)
# ========================================

# Define candidate loci (same as before)
candidates <- data.frame(
  chromosome = c("LG1", "LG2", "LG3"), 
  position = c(12345, 67890, 111213)
)

# OLD WAY (complex):
# go_data <- extract_go_terms_for_enrichment(con, foreground_file_id, background_file_id)
# bp_results <- perform_go_enrichment(go_data, "BP")
# mf_results <- perform_go_enrichment(go_data, "MF")
# # ... many steps

# NEW WAY (simple):

# GO enrichment analysis - all ontologies at once
go_results <- run_go_enrichment_analysis(
  annotations, 
  candidate_loci = candidates,
  ontologies = c("BP", "MF", "CC"),
  significance_threshold = 0.05
)

# KEGG pathway enrichment analysis
kegg_results <- run_kegg_enrichment_analysis(
  annotations,
  candidate_loci = candidates,
  significance_threshold = 0.05
)

# ========================================
# STEP 4: Explore Results
# ========================================

# View GO enrichment results
print(go_results$BP)        # Biological Process results
print(go_results$MF)        # Molecular Function results  
print(go_results$CC)        # Cellular Component results

# View KEGG pathway results
print(kegg_results)

# Count significant results
sapply(go_results, function(x) sum(x$p_adjusted < 0.05, na.rm = TRUE))
sum(kegg_results$p_adjusted < 0.05, na.rm = TRUE)

# ========================================
# STEP 5: Export and Share
# ========================================

# The annotations data frame can be used for:

# 1. Custom analysis in R
pathway_counts <- table(unlist(strsplit(annotations$kegg_modules, ";")))
print(sort(pathway_counts, decreasing = TRUE))

# 2. Export to Excel/CSV for collaborators
write.csv(annotations, "complete_functional_annotations.csv", row.names = FALSE)

# 3. Filter for specific genes or pathways
metabolism_loci <- annotations[grepl("Metabolism", annotations$kegg_modules), ]

# 4. Integration with other tools (already in standard format)
# Can be directly imported into Cytoscape, IPA, DAVID, etc.

# ========================================
# Key Benefits of New Workflow
# ========================================

# ✅ SIMPLIFIED: 3 main functions instead of 10+
# ✅ STANDARDIZED: One consistent data format  
# ✅ ENHANCED: KEGGREST and BRITE always applied
# ✅ EXPORTABLE: CSV output for external tools
# ✅ POWERFUL: All statistical capabilities preserved
# ✅ FUTURE-PROOF: Easy to add new annotation types

# Clean up
close_funseq_db(con)