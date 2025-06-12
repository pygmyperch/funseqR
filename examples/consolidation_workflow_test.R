# Complete Workflow Test for Scaffold Consolidation
# This script demonstrates the proper setup and usage of scaffold consolidation

library(funseqR)

cat("=== Scaffold Consolidation Workflow Test ===\n\n")

# For this test, you need to:
# 1. First define main chromosomes in your database  
# 2. Then the summary functions will use consolidated view by default

cat("Step 1: Define main chromosomes (run this once per database)\n")
cat("-------------------------------------------------------\n")
cat("# Connect to your database\n")
cat("con <- connect_funseq_db('your_database.db')\n\n")

cat("# Define which chromosomes to keep separate\n")
cat("main_chroms <- c('LG1', 'LG2', 'LG3', 'LG4', 'LG5', 'LG6', 'LG7', 'LG8',\n")
cat("                 'LG9', 'LG10', 'LG11', 'LG12', 'LG13', 'LG14', 'LG15',\n") 
cat("                 'LG16', 'LG17', 'LG18', 'LG19', 'LG20', 'LG21', 'LG22',\n")
cat("                 'LG23', 'LG24')\n\n")

cat("# Store this configuration in the database\n")
cat("define_chromosomes(con, main_chroms)\n")
cat("# Output: Main chromosomes defined: LG1, LG2, LG3, LG4, LG5, ...\n")
cat("#         All other scaffolds will be grouped as 'US' in visualizations\n\n")

cat("Step 2: Use consolidated visualizations (automatic after step 1)\n")
cat("--------------------------------------------------------------\n")
cat("# Now summary functions will automatically use consolidated view\n")
cat("get_database_summary(con)  # Shows LG1-LG24 + US instead of 300+ scaffolds\n\n")

cat("# Or override to see all original scaffolds\n")
cat("get_database_summary(con, use_consolidated_chrom = FALSE)  # Shows all scaffolds\n\n")

cat("Step 3: Manual consolidation for custom analysis\n")
cat("-----------------------------------------------\n")
cat("# Get VCF data and manually consolidate\n")
cat("vcf_data <- get_vcf_data(con, file_id = 1)\n")
cat("consolidated <- consolidate_scaffolds(vcf_data, main_chroms)\n")
cat("table(consolidated$chromosome)  # Clean summary\n\n")

cat("=== To Fix Your Current Issue ===\n")
cat("You need to run this command on your actual database:\n\n")
cat("main_chroms <- c('LG1', 'LG2', 'LG3', 'LG4', 'LG5', 'LG6', 'LG7', 'LG8',\n")
cat("                 'LG9', 'LG10', 'LG11', 'LG12', 'LG13', 'LG14', 'LG15',\n")
cat("                 'LG16', 'LG17', 'LG18', 'LG19', 'LG20', 'LG21', 'LG22',\n")
cat("                 'LG23', 'LG24')\n")
cat("define_chromosomes(con, main_chroms)\n\n")

cat("After running this, get_database_summary(con) should show:\n")
cat("- Consolidated view with LG1-LG24 + US in the printed output\n")
cat("- chromosome_distribution: original data (312 chromosomes)\n") 
cat("- chromosome_distribution_display: consolidated data (25 chromosomes)\n\n")

cat("The key insight: consolidation only works AFTER you define main chromosomes!\n")