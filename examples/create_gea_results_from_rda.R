#' Create GEA Results DataFrame from RDA Q-values
#'
#' This script demonstrates how to create a gea_results dataframe for Manhattan plots
#' by combining RDA q-values with genomic coordinates from the database.

library(funseqR)

# Connect to database
con <- connect_funseq_db("funseq_project.db")

# Your VCF file and import (already done)
vcf_file <- "SA448_14699.vcf"
# vcf_import <- import_vcf_to_db(con, vcf_file)

# Get the file_id for your VCF file
vcf_files <- DBI::dbGetQuery(con, "
  SELECT file_id, file_name
  FROM input_files
  WHERE file_type = 'vcf' AND file_name LIKE '%SA448_14699%'
")

if (nrow(vcf_files) == 0) {
  stop("VCF file SA448_14699.vcf not found in database. Please import it first.")
}

file_id <- vcf_files$file_id[1]
cat("Using VCF file ID:", file_id, "\n")

# Get genomic coordinates from database in the same order as VCF
vcf_coords <- DBI::dbGetQuery(con, "
  SELECT vcf_id, chromosome, position, ref, alt
  FROM vcf_data
  WHERE file_id = ?
  ORDER BY vcf_id
", list(file_id))

cat("Retrieved", nrow(vcf_coords), "variants from database\n")

# Assuming you have your RDA results loaded
# rda.simple.pq <- your_rda_results  # Replace with your actual RDA data

# Check that the lengths match
if (exists("rda.simple.pq")) {
  if (length(rda.simple.pq$q.values) != nrow(vcf_coords)) {
    stop("Length mismatch: RDA has ", length(rda.simple.pq$q.values),
         " q-values but VCF has ", nrow(vcf_coords), " variants")
  }

  # Create gea_results dataframe
  gea_results <- data.frame(
    vcf_id = vcf_coords$vcf_id,
    chromosome = vcf_coords$chromosome,
    position = vcf_coords$position,
    ref = vcf_coords$ref,
    alt = vcf_coords$alt,
    q_value = rda.simple.pq$q.values,
    p_value = rda.simple.pq$q.values,  # Using q-values as p-values for now
    stringsAsFactors = FALSE
  )

  # Add -log10 transformed values for Manhattan plots
  gea_results$neg_log10_q <- -log10(gea_results$q_value)
  gea_results$neg_log10_p <- -log10(gea_results$p_value)

  # Add significance thresholds
  gea_results$significant_005 <- gea_results$q_value < 0.05
  gea_results$significant_001 <- gea_results$q_value < 0.01

  cat("Created gea_results with", nrow(gea_results), "variants\n")
  cat("Significant variants (q < 0.05):", sum(gea_results$significant_005), "\n")
  cat("Highly significant variants (q < 0.01):", sum(gea_results$significant_001), "\n")

  # Preview the results
  cat("\nFirst few rows of gea_results:\n")
  print(head(gea_results))

  # Save for later use
  write.table(gea_results, "gea_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  cat("\nSaved gea_results to 'gea_results.txt'\n")

} else {
  cat("Please load your rda.simple.pq data first:\n")
  cat("Example:\n")
  cat("# Load your RDA results\n")
  cat("# rda.simple.pq <- readRDS('path/to/your/rda_results.rds')\n")
  cat("# or load('path/to/your/rda_results.RData')\n")
  cat("\n")

  # Create example structure for reference
  cat("Your rda.simple.pq should have a structure like:\n")
  cat("rda.simple.pq$q.values with", nrow(vcf_coords), "values\n")

  # Create example gea_results structure
  example_gea_results <- data.frame(
    vcf_id = vcf_coords$vcf_id,
    chromosome = vcf_coords$chromosome,
    position = vcf_coords$position,
    ref = vcf_coords$ref,
    alt = vcf_coords$alt,
    q_value = runif(nrow(vcf_coords), 0.001, 0.5),  # Example random q-values
    p_value = runif(nrow(vcf_coords), 0.001, 0.5),  # Example random p-values
    stringsAsFactors = FALSE
  )

  example_gea_results$neg_log10_q <- -log10(example_gea_results$q_value)
  example_gea_results$neg_log10_p <- -log10(example_gea_results$p_value)
  example_gea_results$significant_005 <- example_gea_results$q_value < 0.05
  example_gea_results$significant_001 <- example_gea_results$q_value < 0.01

  cat("\nExample gea_results structure:\n")
  print(str(example_gea_results))

  # Assign to global environment for user reference
  assign("gea_results", example_gea_results, envir = .GlobalEnv)
  cat("\nCreated example 'gea_results' with random values for testing\n")
}

# Alternative approach if you have p-values instead of q-values
create_gea_from_pvalues <- function(vcf_coords, p_values) {
  if (length(p_values) != nrow(vcf_coords)) {
    stop("Length mismatch between p-values and VCF coordinates")
  }

  # Calculate q-values using FDR correction
  q_values <- p.adjust(p_values, method = "fdr")

  gea_results <- data.frame(
    vcf_id = vcf_coords$vcf_id,
    chromosome = vcf_coords$chromosome,
    position = vcf_coords$position,
    ref = vcf_coords$ref,
    alt = vcf_coords$alt,
    p_value = p_values,
    q_value = q_values,
    neg_log10_p = -log10(p_values),
    neg_log10_q = -log10(q_values),
    significant_005 = q_values < 0.05,
    significant_001 = q_values < 0.01,
    stringsAsFactors = FALSE
  )

  return(gea_results)
}

# Example usage:
# If you have p-values: gea_results <- create_gea_from_pvalues(vcf_coords, your_p_values)

cat("\nTo use this data with Manhattan plots:\n")
cat("# Example Manhattan plot call (once you have functional Manhattan plot function):\n")
cat("# manhattan_plot <- create_functional_manhattan_plot(\n")
cat("#   con, enrichment_results, candidate_file_id,\n")
cat("#   gea_results, y_metric = 'neg_log10_q',\n")
cat("#   title = 'RDA Q-values Manhattan Plot'\n")
cat("# )\n")

close_funseq_db(con)
