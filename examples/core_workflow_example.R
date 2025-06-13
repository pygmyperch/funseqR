# funseqR Core Workflow Example
# 
# This example demonstrates the complete core workflow for functional
# annotation of genomic variants using funseqR
#
# Workflow:
# 1. Create/connect to database
# 2. Import VCF and reference genome data
# 3. Extract flanking sequences around variants
# 4. Perform BLAST searches against protein databases
# 5. Annotate BLAST results with UniProt data
# 6. Perform GO enrichment analysis (optional)

library(funseqR)

# Define paths - Update these for your data
db_path <- "funseq_project.db"
vcf_file <- "your_variants.vcf"
vcf_file_cand <- "your_candidate_variants.vcf"  # For GO enrichment
ref_genome_file <- "your_reference_genome.fna"
blast_db_path <- "/path/to/your/blast/database"
blast_db_name <- "your_db_name"

# Define main chromosomes (update for your species)
main_chroms <- c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", 
                 "LG9", "LG10", "LG11", "LG12", "LG13", "LG14", "LG15", 
                 "LG16", "LG17", "LG18", "LG19", "LG20", "LG21", "LG22", 
                 "LG23", "LG24")

# === STEP 1: Create/Connect Database ===
cat("Creating funseqR database...\n")
con <- create_funseq_db(db_path, force = TRUE)

# Alternatively, connect to existing database:
# con <- connect_funseq_db(db_path, verbose = TRUE)

# === STEP 2: Import Data ===
cat("Importing VCF data...\n")
vcf_import <- import_vcf_to_db(con, vcf_file)

# Import candidate variants (optional, for GO enrichment)
# vcf_cand_import <- import_vcf_to_db(con, vcf_file_cand)

cat("Importing reference genome...\n")
ref_import <- import_reference_to_db(con, ref_genome_file, 
                                     genome_name = "Your_Species", 
                                     genome_build = "v1.0")

# Define main chromosomes for consolidation
cat("Defining main chromosomes...\n")
define_chromosomes(con, main_chroms)

# === STEP 3: Extract Flanking Sequences ===
cat("Extracting flanking sequences...\n")
flanking_import <- import_flanking_seqs_to_db(
  con,
  vcf_import$file_id,
  ref_import$genome_id,
  flank_size = 300,
  translate_flanks = TRUE,  # Generate ORF translations
  threads = 4,
  batch_size = 1000,
  verbose = TRUE
)

# Check database status
get_database_summary(con)

# === STEP 4: Perform BLAST Searches ===
cat("Performing BLAST search...\n")

# BLAST ORF translations against protein database
blast_results <- perform_blast_db(
  con,
  vcf_import$file_id,
  db_path = blast_db_path,
  db_name = blast_db_name,
  blast_type = "diamond_blastx",  # Fast DIAMOND BLAST
  seq_type = "orf_nuc",          # Use ORF nucleotide sequences
  e_value = 1e-10,
  max_hits = 5,
  threads = 4
)

# Optional: BLAST raw sequences as well
# blast_results_raw <- perform_blast_db(
#   con,
#   vcf_import$file_id,
#   db_path = blast_db_path,
#   db_name = blast_db_name,
#   blast_type = "diamond_blastx",
#   seq_type = "raw",
#   e_value = 1e-10,
#   max_hits = 5,
#   threads = 4
# )

# Check BLAST results
cat("BLAST hits found:", count_blast_results(con, blast_results$blast_param_id), "\n")

# === STEP 5: Annotate BLAST Results ===
cat("Annotating BLAST results with UniProt data...\n")
annotation_results <- annotate_blast_results(
  con,
  blast_results$blast_param_id,
  evidence_keep = c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "TAS", "IC", "IEA"),
  verbose = TRUE
)

# Check annotation results
cat("Annotations created:", annotation_results$annotated_results, "\n")
cat("GO terms extracted:", annotation_results$go_terms, "\n")

# Get detailed annotations
annotations <- get_annotations(con, blast_results$blast_param_id)
cat("Retrieved", nrow(annotations$annotations), "annotations\n")
cat("Retrieved", nrow(annotations$go_terms), "GO terms\n")

# === STEP 6: GO Enrichment Analysis (Optional) ===
# Requires candidate variants file

# if (exists("vcf_cand_import")) {
#   cat("Performing GO enrichment analysis...\n")
#   
#   # Run complete GO enrichment workflow
#   enrichment_results <- run_go_enrichment_workflow(
#     con,
#     candidate_vcf_file = vcf_file_cand,
#     background_file_id = vcf_import$file_id,
#     ontologies = c("BP", "MF"),
#     create_plots = TRUE,
#     store_results = TRUE,
#     verbose = TRUE
#   )
#   
#   # Check results
#   if (enrichment_results$status == "success") {
#     cat("GO enrichment completed successfully!\n")
#     print(enrichment_results$summary$summary_stats)
#   }
# }

# === STEP 7: Database Summary ===
cat("Final database summary:\n")
final_summary <- get_database_summary(con, verbose = TRUE)

# Close database connection
close_funseq_db(con)

cat("Core workflow completed successfully!\n")

# === Expected Output ===
# Your database should now contain:
# - VCF variant data
# - Reference genome sequences  
# - Flanking sequences around variants (both raw and ORF-translated)
# - BLAST search results
# - UniProt functional annotations
# - GO terms associated with annotations
# - (Optional) GO enrichment results comparing candidate vs background variants