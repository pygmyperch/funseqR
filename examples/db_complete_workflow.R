# Complete workflow using funseqR with SQLite database integration
#
# This script demonstrates a full end-to-end workflow using the database-integrated
# functions for improved reproducibility and data management.

library(funseqR)
library(Biostrings)
library(ggplot2)

HOME <- "/Users/brau0037/Library/CloudStorage/GoogleDrive-pygmyperch@gmail.com/My Drive/git_repos/temp"
setwd(HOME)


# Define paths
db_path <- "funseq_project.db"
vcf_file <- "SA448_14699.vcf"
ref_genome_file <- "Chrysophrys_auratus.v.1.0.all.assembly.units.fna"
blast_db_path <- "/Volumes/SSD2TB/blastDBs/teleost/"
blast_db_name <- "teleostei_db"

# Step 1: Create a new database
con <- create_funseq_db(db_path, force = TRUE)
cat("Created new funseqR database:", db_path, "\n")

# Step 2: Create a new project
project_id <- create_project(
  con,
  project_name = "SnapperAnalysis",
  description = "Analysis of snapper SNPs with functional annotations"
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

# Step 7: Process BLAST results and annotate
annotation_result <- annotate_blast_results(
  con,
  blast_result$blast_param_id,
  max_hits = 3,
  e_value_threshold = 1e-10,
  delay = 1
)
cat("Annotation summary:\n")
cat("- Unique accessions:", annotation_result$unique_accessions, "\n")
cat("- Successful extractions:", annotation_result$successful_extractions, "\n")
cat("- Annotated results:", annotation_result$annotated_results, "\n")
cat("- GO terms:", annotation_result$go_terms, "\n")
cat("- KEGG references:", annotation_result$kegg_refs, "\n")

# Step 8: Retrieve and analyze annotations
annotations <- get_annotations(
  con,
  blast_result$blast_param_id,
  include_go = TRUE,
  include_kegg = TRUE,
  include_vcf_info = TRUE
)

# Display annotation statistics
cat("\nAnnotation Statistics:\n")
cat("Total annotations:", nrow(annotations$annotations), "\n")
cat("Total GO terms:", nrow(annotations$go_terms), "\n")
cat("Total KEGG references:", nrow(annotations$kegg_refs), "\n")

# Analyze GO terms by category
if (nrow(annotations$go_terms) > 0) {
  go_categories <- table(annotations$go_terms$go_category)
  cat("\nGO terms by category:\n")
  print(go_categories)

  # Create a directory for outputs
  output_dir <- "funseq_output"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  # Create a bar plot of GO categories
  if (require(ggplot2)) {
    go_cat_df <- data.frame(
      Category = c("Biological Process", "Molecular Function", "Cellular Component"),
      Code = c("P", "F", "C"),
      Count = c(
        sum(annotations$go_terms$go_category == "P"),
        sum(annotations$go_terms$go_category == "F"),
        sum(annotations$go_terms$go_category == "C")
      )
    )

    p <- ggplot(go_cat_df, aes(x = Category, y = Count, fill = Category)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      labs(
        title = "GO Terms by Category",
        x = "GO Category",
        y = "Number of Terms"
      )

    # Save the plot
    ggsave(file.path(output_dir, "go_categories.png"), p, width = 8, height = 6)
    cat("Created GO categories plot in", file.path(output_dir, "go_categories.png"), "\n")
  }

  # Top 10 most frequent GO terms
  top_go_terms <- as.data.frame(table(annotations$go_terms$go_term))
  colnames(top_go_terms) <- c("GO_Term", "Frequency")
  top_go_terms <- top_go_terms[order(-top_go_terms$Frequency), ]
  top_go_terms <- head(top_go_terms, 10)

  cat("\nTop 10 most frequent GO terms:\n")
  print(top_go_terms)

  # Save GO terms to a file
  write.csv(
    annotations$go_terms[, c("go_id", "go_term", "go_category", "go_evidence", "uniprot_accession")],
    file = file.path(output_dir, "go_terms.csv"),
    row.names = FALSE
  )
  cat("Saved GO terms to", file.path(output_dir, "go_terms.csv"), "\n")
}

# Analyze annotations by chromosome
if (nrow(annotations$annotations) > 0) {
  annotations_by_chrom <- table(annotations$annotations$chromosome)
  cat("\nAnnotations by chromosome:\n")
  print(annotations_by_chrom)

  # Export annotation data
  write.csv(
    annotations$annotations[, c("chromosome", "position", "uniprot_accession", "gene_names", "e_value", "bit_score")],
    file = file.path(output_dir, "annotations.csv"),
    row.names = FALSE
  )
  cat("Saved annotations to", file.path(output_dir, "annotations.csv"), "\n")

  # Export a summary file linking SNPs to genes
  snp_gene_summary <- annotations$annotations[!duplicated(annotations$annotations$vcf_id), ]
  snp_gene_summary <- snp_gene_summary[, c("chromosome", "position", "ref", "alt", "uniprot_accession", "gene_names", "e_value")]
  write.csv(
    snp_gene_summary,
    file = file.path(output_dir, "snp_gene_summary.csv"),
    row.names = FALSE
  )
  cat("Saved SNP-gene summary to", file.path(output_dir, "snp_gene_summary.csv"), "\n")
}

# Step 9: Export a complete data report
if (require(rmarkdown) && file.exists("annotation_report_template.Rmd")) {
  # Generate a report using an R Markdown template
  rmarkdown::render(
    "annotation_report_template.Rmd",
    output_file = file.path(output_dir, "annotation_report.html"),
    params = list(
      project_id = project_id,
      blast_param_id = blast_result$blast_param_id,
      vcf_file_id = vcf_import$file_id,
      db_path = db_path
    )
  )
  cat("Generated annotation report at", file.path(output_dir, "annotation_report.html"), "\n")
}

# Step 10: Record workflow in database for reproducibility
workflow_description <- paste(
  "Full workflow: VCF import -> Reference genome -> Flanking sequences -> BLAST -> Annotation",
  "\nVCF file:", vcf_file,
  "\nReference genome:", ref_genome_file,
  "\nBLAST database:", file.path(blast_db_path, blast_db_name),
  "\nBLAST type: blastx, E-value:", 1e-10, ", Max hits:", 5,
  "\nAnnotation: Max hits per query:", 3, ", E-value threshold:", 1e-10
)

DBI::dbExecute(
  con,
  "INSERT INTO metadata (key, value) VALUES (?, ?)",
  params = list(
    paste0("workflow_", format(Sys.time(), "%Y%m%d_%H%M%S")),
    workflow_description
  )
)

# Close the database connection
close_funseq_db(con)
cat("\nDatabase connection closed\n")
cat("Workflow completed successfully!\n")



# Load the package
library(funseqR)
library(Biostrings)
library(ggplot2)

HOME <- "/Users/brau0037/Library/CloudStorage/GoogleDrive-pygmyperch@gmail.com/My Drive/git_repos/temp"
setwd(HOME)

# Reconnect to the existing database
# Replace "funseq_example.db" with your actual database path
con <- connect_funseq_db("funseq_project.db", verbose = TRUE)

# List available projects to find your project_id
projects <- list_projects(con)
print(projects)
project_id <- projects$project_id[1]  # Adjust as needed to select your project

# List VCF files in the project
vcf_files <- list_input_files(con, project_id)
print(vcf_files[vcf_files$file_type == "vcf", ])
vcf_file_id <- vcf_files$file_id[vcf_files$file_type == "vcf"][1]  # Adjust as needed

# List BLAST parameter sets if you need to continue with annotation
blast_params <- list_blast_params(con, project_id)
print(blast_params)
blast_param_id <- blast_params$blast_param_id[1]  # Adjust as needed

# Now you can continue with your workflow, e.g., annotation
annotation_result <- annotate_blast_results(
  con,
  blast_param_id,
  max_hits = 1,
  e_value_threshold = 1e-10,
  batch_size = 400,
  delay = 3
)

# When done, close the connection
close_funseq_db(con)
