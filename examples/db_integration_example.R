# Example workflow using funseqR with SQLite database integration
#
# This script demonstrates how to use the database-integrated functions
# to store and manage data in a funseqR workflow.

library(funseqR)
library(Biostrings)

# Define paths
db_path <- "funseq_example.db"
vcf_file <- "SA448_14699.vcf"
ref_genome_file <- "Chrysophrys_auratus.v.1.0.all.assembly.units.fna"

# Step 1: Create a new database
con <- create_funseq_db(db_path, force = TRUE)
cat("Created new funseqR database:", db_path, "\n")

# Step 2: Create a new project
project_id <- create_project(
  con,
  project_name = "SnapperAnalysis",
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

# Step 6: Get VCF data from database
vcf_data <- get_vcf_data(con, vcf_import$file_id, limit = 10)
cat("First 10 VCF entries:\n")
print(vcf_data[, c("chromosome", "position", "ref", "alt")])

# Step 7: Get flanking sequences
flanking_seqs <- get_flanking_sequences(con, vcf_import$file_id, as_dna_string_set = TRUE)
cat("Retrieved", length(flanking_seqs), "flanking sequences\n")
cat("First 3 flanking sequence names:\n")
print(names(flanking_seqs)[1:3])

# Step 8: Convert VCF to BED format
bed_data <- vcf2bed_db(con, vcf_import$file_id)
cat("First 5 BED entries:\n")
print(head(bed_data, 5))

# Step 9: List all input files for the project
files <- list_input_files(con, project_id)
cat("Files in project:\n")
print(files[, c("file_id", "file_type", "file_name")])

# Step A: Use integrated database information for BLAST search
# (This is a non-functional example, but shows how the workflow would continue)
#
# Instead of this:
# blast_results <- perform_blast(
#   flanking_seqs,
#   db_path = "/path/to/blast/db",
#   db_name = "teleostei_db",
#   blast_type = "blastx",
#   e_value = 1e-10,
#   max_hits = 5,
#   threads = 4,
#   output_base = "snapper_all_SwPr_teleost"
# )
#
# You would eventually do this:
# blast_results <- perform_blast_db(
#   con,
#   vcf_file_id = vcf_import$file_id,
#   db_path = "/path/to/blast/db",
#   db_name = "teleostei_db",
#   blast_type = "blastx",
#   e_value = 1e-10,
#   max_hits = 5,
#   threads = 4
# )

# Close the database connection
close_funseq_db(con)
cat("Database connection closed\n")
