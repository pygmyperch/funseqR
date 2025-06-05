# Test BLAST on just LG2 sequences

library(funseqR)
library(DBI)
library(Biostrings)

# Function to run BLAST on specific chromosome
test_chromosome_blast <- function(con, chromosome, blast_db_path, blast_db_name, 
                                 project_id = 1, vcf_file_id = 1) {
  
  # Get flanking sequences for specific chromosome
  cat("Getting flanking sequences for", chromosome, "...\n")
  
  flanking_query <- "
    SELECT f.flanking_id, f.sequence, v.chromosome, v.position
    FROM flanking_sequences f
    JOIN vcf_data v ON f.vcf_id = v.vcf_id
    WHERE v.chromosome = ? AND v.file_id = ?
    ORDER BY v.position
  "
  
  flanking_data <- DBI::dbGetQuery(
    con, 
    flanking_query, 
    params = list(chromosome, vcf_file_id)
  )
  
  if (nrow(flanking_data) == 0) {
    cat("No flanking sequences found for", chromosome, "\n")
    return(NULL)
  }
  
  cat("Found", nrow(flanking_data), "flanking sequences for", chromosome, "\n")
  
  # Convert to DNAStringSet
  sequences <- Biostrings::DNAStringSet(flanking_data$sequence)
  names(sequences) <- paste0("seq_", seq_along(sequences))
  
  # Create output files
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  output_base <- paste0("test_", chromosome, "_", timestamp)
  output_fasta <- paste0(output_base, "_query_seqs.fasta")
  output_blast <- paste0(output_base, "_blastx.txt")
  
  # Write sequences to FASTA file
  Biostrings::writeXStringSet(sequences, output_fasta)
  cat("Query sequences written to:", output_fasta, "\n")
  
  # Construct BLAST command
  db_file <- file.path(blast_db_path, blast_db_name)
  outfmt <- "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
  
  blast_command <- paste(
    "blastx",
    "-db", shQuote(db_file),
    "-query", shQuote(output_fasta),
    "-out", shQuote(output_blast),
    "-outfmt", shQuote(outfmt),
    "-evalue", "1e-10",
    "-max_target_seqs", "5",
    "-num_threads", "4"
  )
  
  # Run BLAST
  cat("Executing BLAST command...\n")
  cat("Command:", blast_command, "\n")
  
  system_result <- system(blast_command)
  
  if (system_result != 0) {
    cat("BLAST command failed with exit status:", system_result, "\n")
    return(NULL)
  }
  
  cat("BLAST completed successfully\n")
  cat("Results written to:", output_blast, "\n")
  
  # Check results
  if (file.exists(output_blast) && file.size(output_blast) > 0) {
    # Read and display first few lines
    blast_results <- read.table(output_blast, stringsAsFactors = FALSE)
    cat("Found", nrow(blast_results), "BLAST hits\n")
    cat("First few results:\n")
    print(head(blast_results))
    
    return(list(
      chromosome = chromosome,
      sequences_processed = nrow(flanking_data),
      blast_hits = nrow(blast_results),
      output_fasta = output_fasta,
      output_blast = output_blast,
      flanking_data = flanking_data,
      blast_results = blast_results
    ))
  } else {
    cat("No BLAST hits found\n")
    return(list(
      chromosome = chromosome,
      sequences_processed = nrow(flanking_data),
      blast_hits = 0,
      output_fasta = output_fasta,
      output_blast = output_blast,
      flanking_data = flanking_data
    ))
  }
}

# Set up connection and paths (modify as needed)
HOME <- "/Users/brau0037/Library/CloudStorage/GoogleDrive-pygmyperch@gmail.com/My Drive/projects/snapperSG/south_coast/annotation/funseq_annotationMS"
setwd(HOME)

db_path <- "funseq_project.db"
blast_db_path <- "/Volumes/SSD2TB/blastDBs/teleost/"
blast_db_name <- "teleostei_db"

# Connect to database
con <- connect_funseq_db(db_path, verbose = TRUE)

# Test LG2
cat("Testing BLAST on LG2...\n")
lg2_result <- test_chromosome_blast(
  con, 
  chromosome = "LG2", 
  blast_db_path = blast_db_path, 
  blast_db_name = blast_db_name
)

# Display results
if (!is.null(lg2_result)) {
  cat("\n=== LG2 BLAST Test Results ===\n")
  cat("Chromosome:", lg2_result$chromosome, "\n")
  cat("Sequences processed:", lg2_result$sequences_processed, "\n")
  cat("BLAST hits found:", lg2_result$blast_hits, "\n")
  cat("Output files:", lg2_result$output_fasta, "and", lg2_result$output_blast, "\n")
}

# Close connection
close_funseq_db(con)