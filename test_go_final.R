#!/usr/bin/env Rscript

# Final test of GO enrichment implementation
library(funseqR)

project_dir <- "/Users/brau0037/Library/CloudStorage/GoogleDrive-pygmyperch@gmail.com/My Drive/projects/snapperSG/south_coast/annotation/funseq_annotationMS"
db_path <- file.path(project_dir, "funseq_project.db")
candidate_vcf <- file.path(project_dir, "SA448_855.vcf")

con <- connect_funseq_db(db_path)

cat("=== Testing GO Enrichment with Snapper Data ===\n")

results <- run_go_enrichment_workflow(
  con = con,
  project_id = 1,
  candidate_vcf_file = candidate_vcf,
  background_file_id = 1,
  ontologies = c("BP", "MF"),  # Test both ontologies
  min_genes = 3,
  store_results = FALSE,
  create_plots = TRUE,
  verbose = TRUE
)

cat("=== RESULTS SUMMARY ===\n")
print(results$summary)

# Check Biological Process results
if (!is.null(results$enrichment_results$BP)) {
  bp_results <- results$enrichment_results$BP
  sig_count <- sum(bp_results$p_adjusted < 0.05, na.rm = TRUE)
  cat("Biological Process - Significant terms:", sig_count, "\n")
  
  if (sig_count > 0) {
    sig_results <- bp_results[bp_results$p_adjusted < 0.05, ]
    sig_results <- sig_results[order(sig_results$p_adjusted), ]
    top5 <- head(sig_results, 5)
    cat("Top 5 enriched BP terms:\n")
    for (i in 1:nrow(top5)) {
      cat(sprintf("  %d. %s (Fold: %.2f, FDR: %.2e)\n", 
                  i, top5$go_term[i], top5$fold_enrichment[i], top5$p_adjusted[i]))
    }
  }
}

# Check Molecular Function results
if (!is.null(results$enrichment_results$MF)) {
  mf_results <- results$enrichment_results$MF
  sig_count <- sum(mf_results$p_adjusted < 0.05, na.rm = TRUE)
  cat("Molecular Function - Significant terms:", sig_count, "\n")
  
  if (sig_count > 0) {
    sig_results <- mf_results[mf_results$p_adjusted < 0.05, ]
    sig_results <- sig_results[order(sig_results$p_adjusted), ]
    top3 <- head(sig_results, 3)
    cat("Top 3 enriched MF terms:\n")
    for (i in 1:nrow(top3)) {
      cat(sprintf("  %d. %s (Fold: %.2f, FDR: %.2e)\n", 
                  i, top3$go_term[i], top3$fold_enrichment[i], top3$p_adjusted[i]))
    }
  }
}

# Test plots
if (!is.null(results$plots$BP_bubble)) {
  cat("BP bubble plot created successfully\n")
}

if (!is.null(results$plots$MF_bubble)) {
  cat("MF bubble plot created successfully\n")
}

close_funseq_db(con)

cat("=== GO ENRICHMENT TEST COMPLETED SUCCESSFULLY ===\n")