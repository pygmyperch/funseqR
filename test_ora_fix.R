#!/usr/bin/env Rscript

# Test script to verify ORA storage fix
setwd('/Users/brau0037/Library/CloudStorage/GoogleDrive-pygmyperch@gmail.com/My Drive/git_repos/funseqR')

library(funseqR)

con <- connect_funseq_db('test_data/snapper_candidates.db')

cat("Testing ORA function with fixed parameter binding...\n")

ORA_results <- ora(con, significance_threshold = 0.1)

cat("Test completed successfully!\n")
cat("ORA results status:", ORA_results$status, "\n")