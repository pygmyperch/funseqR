# Script to fix the BLAST mapping issue and re-import results correctly

library(funseqR)
library(DBI)

# Load the fixed import function
source("/Users/brau0037/Library/CloudStorage/GoogleDrive-pygmyperch@gmail.com/My Drive/git_repos/funseqR/R/db_blast_fixed.R")

# Set up paths
HOME <- "/Users/brau0037/Library/CloudStorage/GoogleDrive-pygmyperch@gmail.com/My Drive/projects/snapperSG/south_coast/annotation/funseq_annotationMS"
setwd(HOME)

db_path <- "funseq_project.db"
blast_results_file <- "project_1_vcf_1_teleostei_db_blastx_20250605_150434_blastx.txt"

# Connect to database
con <- connect_funseq_db(db_path, verbose = TRUE)

# Check current state
cat("=== Current Database State ===\n")
blast_params <- DBI::dbGetQuery(con, "SELECT * FROM blast_parameters")
print(blast_params)

current_results <- DBI::dbGetQuery(
  con, 
  "SELECT blast_param_id, COUNT(*) as count FROM blast_results GROUP BY blast_param_id"
)
cat("\nCurrent BLAST results count by parameter ID:\n")
print(current_results)

# Show current chromosome distribution for parameter ID 2
cat("\nCurrent annotations by chromosome (parameter ID 2):\n")
chrom_dist <- DBI::dbGetQuery(
  con,
  "SELECT v.chromosome, COUNT(br.blast_result_id) as blast_results 
   FROM vcf_data v 
   LEFT JOIN flanking_sequences f ON v.vcf_id = f.vcf_id 
   LEFT JOIN blast_results br ON f.flanking_id = br.flanking_id 
   WHERE br.blast_param_id = 2 AND v.chromosome LIKE 'LG%'
   GROUP BY v.chromosome 
   ORDER BY v.chromosome"
)
print(chrom_dist)

# Test the fixed import function
cat("\n=== Testing Fixed Import Function ===\n")

# Re-import BLAST results using the fixed mapping logic
cat("Re-importing BLAST results with fixed mapping...\n")
new_count <- reimport_blast_results(
  con, 
  blast_param_id = 2,  # Use existing parameter ID 2
  results_file = blast_results_file,
  vcf_file_id = 1,
  clear_existing = TRUE,  # Clear the incorrectly mapped results first
  verbose = TRUE
)

cat("\n=== Results After Fix ===\n")
cat("Successfully imported", new_count, "BLAST results\n")

# Check new chromosome distribution
cat("\nNew annotations by chromosome (parameter ID 2):\n")
new_chrom_dist <- DBI::dbGetQuery(
  con,
  "SELECT v.chromosome, COUNT(br.blast_result_id) as blast_results 
   FROM vcf_data v 
   LEFT JOIN flanking_sequences f ON v.vcf_id = f.vcf_id 
   LEFT JOIN blast_results br ON f.flanking_id = br.flanking_id 
   WHERE br.blast_param_id = 2 AND v.chromosome LIKE 'LG%'
   GROUP BY v.chromosome 
   ORDER BY v.chromosome"
)
print(new_chrom_dist)

# Count total linkage groups with results
lg_with_results <- sum(new_chrom_dist$blast_results > 0)
cat("\nLinkage groups with BLAST results:", lg_with_results, "out of 25\n")

# Show some example mappings
cat("\nExample BLAST result mappings:\n")
examples <- DBI::dbGetQuery(
  con,
  "SELECT v.chromosome, v.position, br.hit_accession, br.e_value
   FROM blast_results br
   JOIN flanking_sequences f ON br.flanking_id = f.flanking_id
   JOIN vcf_data v ON f.vcf_id = v.vcf_id
   WHERE br.blast_param_id = 2
   ORDER BY v.chromosome, v.position
   LIMIT 10"
)
print(examples)

# Close connection
close_funseq_db(con)

cat("\n=== Fix Complete ===\n")
cat("You can now run get_annotations() with blast_param_id = 2 to get the corrected results\n")