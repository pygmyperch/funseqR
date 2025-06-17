#' Descriptive Analysis Example for funseqR
#'
#' This example demonstrates how to use the comprehensive descriptive functions
#' to analyze functional annotations, coverage statistics, and generate detailed
#' profiles of candidate loci.

library(funseqR)

cat("=== funseqR Descriptive Analysis Example ===\n\n")

# Connect to database
cat("Step 1: Connecting to database...\n")
# con <- connect_funseq_db("funseq_project.db")

# Note: Replace the above line with your actual database connection
# For this example, we'll show the expected usage patterns

cat("=== Example Usage Patterns ===\n\n")

cat("1. GO ANNOTATION ANALYSIS\n")
cat("-------------------------\n")
cat("# Comprehensive GO term analysis\n")
cat("go_summary <- summarize_go_annotations(con)\n")
cat("print(go_summary$overview)  # Summary by GO category (BP, MF, CC)\n")
cat("print(go_summary$top_terms)  # Most frequent GO terms\n")
cat("print(go_summary$evidence_summary)  # GO evidence codes\n\n")

cat("# Focus on specific BLAST parameter set\n")
cat("go_summary_blast1 <- summarize_go_annotations(con, blast_param_id = 1)\n\n")

cat("# Increase minimum frequency threshold for cleaner results\n")
cat("go_summary_filtered <- summarize_go_annotations(con, min_frequency = 5)\n\n")

cat("2. UNIPROT ANNOTATION ANALYSIS\n")
cat("-------------------------------\n")
cat("# Comprehensive UniProt analysis\n")
cat("uniprot_summary <- summarize_uniprot_annotations(con)\n")
cat("print(uniprot_summary$overview)  # Coverage statistics\n")
cat("print(uniprot_summary$top_proteins)  # Most frequent proteins\n")
cat("print(uniprot_summary$gene_name_stats)  # Gene name patterns\n")
cat("print(uniprot_summary$annotation_quality)  # Quality metrics\n\n")

cat("# Get top 50 most frequent proteins\n")
cat("uniprot_summary_extended <- summarize_uniprot_annotations(con, top_n = 50)\n\n")

cat("3. COMPREHENSIVE FUNCTIONAL PROFILING\n")
cat("-------------------------------------\n")
cat("# Profile all annotated loci in the database\n")
cat("full_profile <- get_functional_profile(con)\n")
cat("print(full_profile$loci_summary)  # Overview of annotation coverage\n")
cat("print(full_profile$go_profile$overview)  # GO term distribution\n")
cat("print(full_profile$protein_profile$overview)  # Protein annotation stats\n")
cat("print(full_profile$functional_diversity)  # Diversity metrics\n\n")

cat("# Profile specific candidate loci from your analysis\n")
cat("# Example: loci from RDA analysis or other GEA methods\n")
cat("candidates <- data.frame(\n")
cat("  chromosome = c('LG1', 'LG2', 'LG3', 'LG1'),\n")
cat("  position = c(123456, 234567, 345678, 456789)\n")
cat(")\n")
cat("candidate_profile <- get_functional_profile(con, candidate_loci = candidates)\n\n")

cat("# Include BLAST quality statistics and KEGG pathways\n")
cat("detailed_profile <- get_functional_profile(\n")
cat("  con, \n")
cat("  candidate_loci = candidates,\n")
cat("  include_blast_stats = TRUE,\n")
cat("  include_pathways = TRUE\n")
cat(")\n")
cat("print(detailed_profile$blast_quality)  # BLAST hit quality\n")
cat("print(detailed_profile$pathway_profile)  # KEGG pathway annotations\n\n")

cat("4. ANNOTATION COVERAGE ANALYSIS\n")
cat("--------------------------------\n")
cat("# Comprehensive coverage statistics\n")
cat("coverage <- get_annotation_coverage(con)\n")
cat("print(coverage$overall)  # Overall coverage across annotation types\n")
cat("print(coverage$by_type)  # Breakdown by GO, UniProt, KEGG\n")
cat("print(coverage$by_chromosome)  # Coverage by chromosome\n")
cat("print(coverage$by_blast_param)  # Coverage by BLAST parameter set\n")
cat("print(coverage$quality_metrics)  # Annotation quality indicators\n\n")

cat("# Faster analysis without chromosome breakdown\n")
cat("coverage_simple <- get_annotation_coverage(\n")
cat("  con, \n")
cat("  by_chromosome = FALSE,\n")
cat("  by_blast_param = FALSE\n")
cat(")\n\n")

cat("=== PRACTICAL ANALYSIS WORKFLOWS ===\n\n")

cat("WORKFLOW 1: Quality Assessment\n")
cat("------------------------------\n")
cat("# 1. Check overall database coverage\n")
cat("coverage <- get_annotation_coverage(con)\n")
cat("cat('Annotation coverage:', coverage$overall$annotation_percent, '%\\n')\n")
cat("cat('GO coverage:', coverage$overall$go_coverage_percent, '%\\n')\n\n")

cat("# 2. Identify best BLAST parameter set\n")
cat("best_blast <- coverage$by_blast_param[which.max(coverage$by_blast_param$total_annotations), ]\n")
cat("cat('Best BLAST parameter:', best_blast$blast_param_id, '\\n')\n\n")

cat("# 3. Assess annotation quality\n")
cat("quality <- coverage$quality_metrics\n")
cat("cat('Loci with rich GO annotation:', quality$loci_with_rich_go_annotation, '\\n')\n")
cat("cat('Average GO terms per locus:', quality$avg_go_terms_per_locus, '\\n')\n\n")

cat("WORKFLOW 2: Candidate Loci Analysis\n")
cat("-----------------------------------\n")
cat("# 1. Load your candidate loci (example from RDA analysis)\n")
cat("# candidates <- get_significant_loci_from_rda(rda_results, threshold = 0.01)\n\n")

cat("# 2. Get comprehensive functional profile\n")
cat("profile <- get_functional_profile(con, candidate_loci = candidates)\n\n")

cat("# 3. Summarize key functional categories\n")
cat("go_categories <- profile$go_profile$overview\n")
cat("cat('Biological Process terms:', go_categories[go_categories$go_category == 'BP', 'unique_go_terms'], '\\n')\n")
cat("cat('Molecular Function terms:', go_categories[go_categories$go_category == 'MF', 'unique_go_terms'], '\\n')\n")
cat("cat('Cellular Component terms:', go_categories[go_categories$go_category == 'CC', 'unique_go_terms'], '\\n')\n\n")

cat("# 4. Identify top functional themes\n")
cat("top_bp_terms <- profile$go_profile$top_terms[profile$go_profile$top_terms$go_category == 'BP', ]\n")
cat("cat('Top Biological Process themes:\\n')\n")
cat("for(i in 1:min(5, nrow(top_bp_terms))) {\n")
cat("  cat('  ', top_bp_terms$go_term[i], ' (', top_bp_terms$frequency[i], ' loci)\\n')\n")
cat("}\n\n")

cat("WORKFLOW 3: Comparative Analysis\n")
cat("--------------------------------\n")
cat("# Compare different BLAST parameter sets\n")
cat("blast1_profile <- get_functional_profile(con, blast_param_id = 1)\n")
cat("blast2_profile <- get_functional_profile(con, blast_param_id = 2)\n\n")

cat("# Compare GO term diversity\n")
cat("cat('BLAST 1 - GO terms:', blast1_profile$functional_diversity$unique_go_terms, '\\n')\n")
cat("cat('BLAST 2 - GO terms:', blast2_profile$functional_diversity$unique_go_terms, '\\n')\n\n")

cat("# Compare protein diversity\n")
cat("cat('BLAST 1 - Proteins:', blast1_profile$functional_diversity$unique_proteins, '\\n')\n")
cat("cat('BLAST 2 - Proteins:', blast2_profile$functional_diversity$unique_proteins, '\\n')\n\n")

cat("=== DATA EXPORT AND REPORTING ===\n\n")

cat("EXPORT SUMMARIES TO FILES\n")
cat("-------------------------\n")
cat("# Export GO summary to CSV\n")
cat("write.csv(go_summary$overview, 'go_summary_overview.csv', row.names = FALSE)\n")
cat("write.csv(go_summary$top_terms, 'go_top_terms.csv', row.names = FALSE)\n\n")

cat("# Export UniProt summary\n")
cat("write.csv(uniprot_summary$overview, 'uniprot_overview.csv', row.names = FALSE)\n")
cat("write.csv(uniprot_summary$top_proteins, 'top_proteins.csv', row.names = FALSE)\n\n")

cat("# Export coverage statistics\n")
cat("write.csv(coverage$overall, 'coverage_overall.csv', row.names = FALSE)\n")
cat("write.csv(coverage$by_type, 'coverage_by_type.csv', row.names = FALSE)\n")
cat("if (!is.null(coverage$by_chromosome)) {\n")
cat("  write.csv(coverage$by_chromosome, 'coverage_by_chromosome.csv', row.names = FALSE)\n")
cat("}\n\n")

cat("# Export functional profile\n")
cat("write.csv(profile$loci_summary, 'candidate_loci_summary.csv', row.names = FALSE)\n")
cat("write.csv(profile$functional_diversity, 'functional_diversity.csv', row.names = FALSE)\n\n")

cat("CREATE SUMMARY REPORT\n")
cat("---------------------\n")
cat("# Generate a comprehensive text report\n")
cat("sink('functional_analysis_report.txt')\n")
cat("cat('=== FUNCTIONAL ANALYSIS REPORT ===\\n\\n')\n")
cat("cat('Generated on:', Sys.Date(), '\\n\\n')\n\n")

cat("cat('DATABASE OVERVIEW\\n')\n")
cat("cat('-----------------\\n')\n")
cat("cat('Total loci:', coverage$overall$total_loci, '\\n')\n")
cat("cat('Annotation coverage:', coverage$overall$annotation_percent, '%\\n')\n")
cat("cat('GO coverage:', coverage$overall$go_coverage_percent, '%\\n')\n\n")

cat("cat('\\nGO ANNOTATION SUMMARY\\n')\n")
cat("cat('---------------------\\n')\n")
cat("print(go_summary$overview)\n\n")

cat("cat('\\nTOP BIOLOGICAL PROCESSES\\n')\n")
cat("cat('------------------------\\n')\n")
cat("bp_terms <- go_summary$top_terms[go_summary$top_terms$go_category == 'BP', ]\n")
cat("print(bp_terms[1:10, ])\n\n")

cat("cat('\\nUNIPROT ANNOTATION QUALITY\\n')\n")
cat("cat('--------------------------\\n')\n")
cat("print(uniprot_summary$annotation_quality)\n\n")

cat("sink()\n")
cat("cat('Report saved to functional_analysis_report.txt\\n')\n\n")

cat("=== ADVANCED ANALYSIS TIPS ===\n\n")

cat("FILTERING AND SUBSETTING\n")
cat("------------------------\n")
cat("# Analyze only high-quality annotations\n")
cat("high_quality_go <- summarize_go_annotations(con, min_frequency = 10)\n\n")

cat("# Focus on specific chromosomes\n")
cat("lg1_candidates <- candidates[candidates$chromosome == 'LG1', ]\n")
cat("lg1_profile <- get_functional_profile(con, candidate_loci = lg1_candidates)\n\n")

cat("# Compare enriched vs background loci\n")
cat("# enriched_loci <- your_enriched_loci_dataframe\n")
cat("# background_loci <- your_background_loci_dataframe\n")
cat("# enriched_profile <- get_functional_profile(con, candidate_loci = enriched_loci)\n")
cat("# background_profile <- get_functional_profile(con, candidate_loci = background_loci)\n\n")

cat("PERFORMANCE OPTIMIZATION\n")
cat("------------------------\n")
cat("# For large datasets, disable detailed breakdowns\n")
cat("fast_coverage <- get_annotation_coverage(\n")
cat("  con, \n")
cat("  by_chromosome = FALSE, \n")
cat("  by_blast_param = FALSE\n")
cat(")\n\n")

cat("# Disable optional analyses for speed\n")
cat("fast_profile <- get_functional_profile(\n")
cat("  con,\n")
cat("  include_blast_stats = FALSE,\n")
cat("  include_pathways = FALSE\n")
cat(")\n\n")

cat("INTEGRATION WITH OTHER ANALYSES\n")
cat("-------------------------------\n")
cat("# Combine with enrichment analysis\n")
cat("# enrichment_results <- run_go_enrichment_workflow(...)\n")
cat("# enriched_loci_summary <- summarize_functional_loci(...)\n")
cat("# enriched_profile <- get_functional_profile(con, candidate_loci = enriched_loci_summary$loci_summary)\n\n")

cat("# Combine with Manhattan plots\n")
cat("# manhattan_plot <- create_functional_manhattan_plot(\n")
cat("#   con, \n")
cat("#   y_values = your_statistical_values,\n")
cat("#   functional_summary = enriched_loci_summary\n")
cat("# )\n\n")

cat("# Close database connection\n")
cat("# close_funseq_db(con)\n\n")

cat("=== SUMMARY OF KEY FUNCTIONS ===\n\n")
cat("1. summarize_go_annotations(con, blast_param_id = NULL, ...)\n")
cat("   - Comprehensive GO term analysis across categories\n")
cat("   - Returns: overview, top_terms, evidence_summary, coverage_stats\n\n")

cat("2. summarize_uniprot_annotations(con, blast_param_id = NULL, ...)\n")
cat("   - UniProt protein annotation analysis\n")
cat("   - Returns: overview, top_proteins, gene_name_stats, annotation_quality\n\n")

cat("3. get_functional_profile(con, candidate_loci = NULL, ...)\n")
cat("   - Comprehensive functional characterization\n")
cat("   - Returns: loci_summary, go_profile, protein_profile, pathway_profile, \n")
cat("             blast_quality, functional_diversity\n\n")

cat("4. get_annotation_coverage(con, by_chromosome = TRUE, ...)\n")
cat("   - Coverage statistics across annotation types and genomic regions\n")
cat("   - Returns: overall, by_type, by_chromosome, by_blast_param, quality_metrics\n\n")

cat("These functions provide the foundation for comprehensive functional\n")
cat("genomics analysis in funseqR. Use them to assess data quality,\n")
cat("characterize candidate loci, and generate detailed reports.\n")