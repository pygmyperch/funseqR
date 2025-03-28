# Example workflow using funseqR with SQLite database integration
# including BLAST search integration
#
# This script demonstrates how to use the database-integrated functions
# to store and manage data in a complete funseqR workflow.

library(funseqR)
library(Biostrings)

# Define paths
db_path <- "funseq_example.db"
vcf_file <- "SA448_14699.vcf"
ref_genome_file <- "Chrysophrys_auratus.v.1.0.all.assembly.units.fna"
blast_db_path <- "/path/to/blastDBs/teleost"
blast_db_name <- "teleostei_db"

# Step 1: Create a new database
con <- create_funseq_db(db_path, force = TRUE)
cat("Created new funseqR database:", db_path, "\n")

# Step 2: Create a new project
project_id <- create_project(
  con,
  project_name = "Snapper Analysis",
  description = "Analysis of snapper data with SNPs"
)
cat("Created new project with ID:", project_id, "\n")

# Step 3: Import VCF data
vcf_import <- import_vcf_to_db(con, project_id, vcf_file)
cat("Imported", vcf_import$vcf_count, "VCF entries with file ID:", vcf_import$file_id, "\n")

# Step 4: Import reference genome
genome_import <- import_reference_to_db(
  con,
  project_id,
  ref_genome_file,
  genome_name = "Chrysophrys_auratus",
  genome_build = "v1.0"
)
cat("Imported reference genome with ID:", genome_import$genome_id, "\n")
cat("Imported", genome_import$sequence_count, "reference sequences\n")

# Step 5: Extract flanking sequences
flanking_import <- import_flanking_seqs_to_db(
  con,
  vcf_import$file_id,
  genome_import$genome_id,
  flank_size = 500
)
cat("Extracted", flanking_import$flanking_count, "flanking sequences\n")

# Step 6: Perform BLAST search
blast_result <- perform_blast_db(
  con,
  project_id,
  vcf_import$file_id,
  db_path = blast_db_path,
  db_name = blast_db_name,
  blast_type = "blastx",
  e_value = 1e-10,
  max_hits = 5,
  threads = 4
)
cat("BLAST search completed with parameter ID:", blast_result$blast_param_id, "\n")
cat("Imported", blast_result$result_count, "BLAST results\n")
cat("Output files stored with base name:", blast_result$output_base, "\n")

# Step 7: Get BLAST results from database
blast_data <- get_blast_results(
  con, 
  blast_result$blast_param_id, 
  e_value_threshold = 1e-10, 
  max_hits_per_query = 3
)
cat("Retrieved", nrow(blast_data), "BLAST results from database\n")

# Display top hits information
if (nrow(blast_data) > 0) {
  cat("\nTop BLAST hits summary:\n")
  top_hits <- blast_data[!duplicated(blast_data$flanking_id), ]
  top_hits_summary <- top_hits[order(top_hits$e_value), ][1:min(10, nrow(top_hits)), ]
  print(top_hits_summary[, c("chromosome", "position", "hit_accession", "percent_identity", "e_value", "bit_score")])
}

# Count results by chromosome
if (nrow(blast_data) > 0) {
  cat("\nBLAST hits by chromosome:\n")
  hits_by_chrom <- table(blast_data$chromosome)
  print(hits_by_chrom)
}

# Step 8: Register a codebase version in your database for reproducibility
DBI::dbExecute(
  con,
  "INSERT INTO metadata (key, value) VALUES (?, ?)",
  params = list(
    paste0("analysis_", format(Sys.time(), "%Y%m%d_%H%M%S")),
    paste0(
      "Analysis completed with funseqR version ", packageVersion("funseqR"),
      ", R version ", R.version$major, ".", R.version$minor,
      ", on ", Sys.info()["nodename"]
    )
  )
)

# Close the database connection
close_funseq_db(con)
cat("Database connection closed\n")

# At a later time, you can reconnect to this database and all your data will be available:
cat("\nReconnecting to database to demonstrate persistence...\n")
con2 <- connect_funseq_db(db_path)

# List projects
projects <- list_projects(con2)
cat("Projects in database:", nrow(projects), "\n")
print(projects[, c("project_id", "project_name", "creation_date")])

# Get BLAST parameters
blast_params <- list_blast_params(con2, project_id)
cat("\nBLAST parameter sets in project:", nrow(blast_params), "\n")
print(blast_params[, c("blast_param_id", "blast_type", "db_name", "e_value", "execution_date")])

# Count BLAST results
result_count <- count_blast_results(con2, project_id = project_id)
cat("\nTotal BLAST results in project:", result_count, "\n")

# Close the second connection
close_funseq_db(con2)
