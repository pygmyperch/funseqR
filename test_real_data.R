#!/usr/bin/env Rscript

# Test script with actual database
library(funseqR)

cat("Testing ORA function with actual database...\n")

con <- connect_funseq_db('funseq_project.db')

cat("Connected to database successfully\n")

ORA_results <- ora(con, significance_threshold = 0.1)

cat("Test completed successfully!\n")
cat("ORA results status:", ORA_results$status, "\n")