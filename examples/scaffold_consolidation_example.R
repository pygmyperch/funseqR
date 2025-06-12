# Scaffold Consolidation Example
# This script demonstrates how to use the new scaffold consolidation functionality
# for cleaner visualizations while maintaining analysis accuracy.

library(funseqR)

# Example workflow using scaffold consolidation
# ============================================

# Step 1: Connect to database (assuming you have VCF data imported)
# con <- connect_funseq_db("my_analysis.db")

# Step 2: Discover available chromosomes
# all_chromosomes <- define_chromosomes(con)
# print(all_chromosomes)
# Example output: [1] "LG1" "LG2" "LG3" "LG4" "LG5" "scaffold_001" "scaffold_002" ... "scaffold_847"

# Step 3: Define main chromosomes (run this once per project)
# main_chromosomes <- c("LG1", "LG2", "LG3", "LG4", "LG5")
# define_chromosomes(con, main_chromosomes)
# Output: Main chromosomes defined: LG1, LG2, LG3, LG4, LG5
#         All other scaffolds will be grouped as 'US' in visualizations

# Step 4: Use consolidated visualizations (default behavior)
# get_database_summary(con)  # Shows: LG1, LG2, LG3, LG4, LG5, US
# get_database_summary(con, use_consolidated_chrom = FALSE)  # Shows all original names

# Step 5: Manual consolidation example
# If you want to consolidate data manually for custom analysis:
# vcf_data <- get_vcf_data(con, file_id = 1)
# consolidated_data <- consolidate_scaffolds(vcf_data, main_chromosomes)
# table(consolidated_data$chromosome)  # Clean summary

# ============================================
# Example with simulated data
# ============================================

# Create example VCF data with many scaffolds
example_vcf <- data.frame(
  vcf_id = 1:1000,
  chromosome = c(
    rep("LG1", 300),
    rep("LG2", 250),
    rep("LG3", 200),
    rep("LG4", 100),
    rep("LG5", 50),
    paste0("scaffold_", sprintf("%03d", 1:100))  # 100 small scaffolds, 1 variant each
  ),
  position = sample(1:1000000, 1000),
  stringsAsFactors = FALSE
)

# Show original distribution (messy)
cat("Original chromosome distribution:\n")
orig_dist <- table(example_vcf$chromosome)
cat("Number of chromosomes:", length(orig_dist), "\n")
cat("Main chromosomes:", sum(orig_dist[c("LG1", "LG2", "LG3", "LG4", "LG5")]), "variants\n")
cat("Small scaffolds:", sum(orig_dist[!names(orig_dist) %in% c("LG1", "LG2", "LG3", "LG4", "LG5")]), "variants\n")

# Apply consolidation
main_chroms <- c("LG1", "LG2", "LG3", "LG4", "LG5")
consolidated_vcf <- consolidate_scaffolds(example_vcf, main_chroms, verbose = TRUE)

# Show clean distribution
cat("\nConsolidated chromosome distribution:\n")
final_dist <- table(consolidated_vcf$chromosome)
print(final_dist)

cat("\nBenefits of consolidation:\n")
cat("- Reduced from", length(orig_dist), "to", length(final_dist), "categories\n")
cat("- Main chromosomes preserved for analysis\n")
cat("- Small scaffolds grouped for cleaner visualization\n")
cat("- Original data unchanged (consolidation only for display)\n")

# ============================================
# Usage recommendations
# ============================================

cat("\n=== Usage Recommendations ===\n")
cat("1. Use define_chromosomes(con) to explore available chromosomes\n")
cat("2. Use define_chromosomes(con, main_chroms) to set project defaults\n")
cat("3. Visualization functions default to consolidated view (use_consolidated_chrom = TRUE)\n")
cat("4. Analysis functions always use original chromosome names\n")
cat("5. Override with use_consolidated_chrom = FALSE to see all scaffolds when needed\n")
cat("6. Manual consolidation available with consolidate_scaffolds() function\n")

# close_funseq_db(con)