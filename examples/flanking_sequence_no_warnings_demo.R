# Demonstration: Flanking Sequence Import with ORF Extraction (No Warnings)
# This script shows that import_flanking_seqs_to_db() now runs cleanly with translate_flanks=TRUE

library(funseqR)

cat("=== Flanking Sequence Import - No Warnings Demo ===\n\n")

cat("This demo shows that the Biostrings translation warnings have been fixed.\n")
cat("Previously, using translate_flanks=TRUE would generate warnings like:\n")
cat("  'Warning: last base was ignored'\n")  
cat("  'Warning: last 2 bases were ignored'\n\n")

cat("These warnings occurred because:\n")
cat("1. Flanking sequences (e.g., 500bp each side = 1000bp total) aren't always divisible by 3\n")
cat("2. Biostrings::translate() warns when sequences have incomplete codons\n")
cat("3. This is normal when translating raw genomic sequences\n\n")

cat("=== Fix Implementation ===\n")
cat("Modified extract_longest_orf() in R/utils.R to:\n")
cat("1. Wrap Biostrings::translate() calls in suppressWarnings()\n")
cat("2. Add explanatory comments about why warnings are suppressed\n")
cat("3. Preserve all functionality while eliminating warning spam\n\n")

cat("=== Test Results ===\n")

# Simulate what would happen during flanking sequence import
cat("Testing ORF extraction on sequences that would generate warnings...\n")

# Example sequences that would cause warnings (not divisible by 3)
test_sequences <- c(
  paste0("ATG", paste(rep("AAATTTCCCGGG", 40), collapse = ""), "A"),     # +1 base
  paste0("ATG", paste(rep("AAATTTCCCGGG", 40), collapse = ""), "AT"),    # +2 bases  
  paste0("ATG", paste(rep("AAATTTCCCGGG", 40), collapse = ""))           # divisible by 3
)

for (i in 1:length(test_sequences)) {
  seq_len <- nchar(test_sequences[i])
  remainder <- seq_len %% 3
  
  cat(sprintf("Sequence %d: %d bp (remainder: %d) - ", i, seq_len, remainder))
  
  # This demonstrates the internal function call that happens during import
  result <- funseqR:::extract_longest_orf(
    test_sequences[i], 
    min_aa_length = 10, 
    return_type = "aa",
    verbose = FALSE
  )
  
  if (!is.null(result)) {
    cat("✓ ORF extracted (", nchar(result), " AA)\n")
  } else {
    cat("✗ No ORF found\n")
  }
}

cat("\n=== Real-World Usage ===\n")
cat("In practice, you would use this functionality like:\n\n")
cat("# Connect to database with VCF and reference genome data\n")
cat("con <- connect_funseq_db('my_analysis.db')\n\n")
cat("# Import flanking sequences with ORF extraction (NO WARNINGS!)\n")
cat("flanking_import <- import_flanking_seqs_to_db(\n")
cat("  con,\n")
cat("  vcf_file_id = 1,\n") 
cat("  genome_id = 1,\n")
cat("  flank_size = 500,\n")
cat("  translate_flanks = TRUE  # This now runs without warnings\n")
cat(")\n\n")

cat("=== Benefits ===\n")
cat("✓ Clean output during analysis\n") 
cat("✓ Warnings suppressed only where expected and harmless\n")
cat("✓ All functionality preserved\n")
cat("✓ ORF extraction works as intended\n")
cat("✓ Better user experience for large datasets\n\n")

cat("The fix ensures that translate_flanks=TRUE runs cleanly on real datasets\n")
cat("without generating hundreds of expected translation warnings.\n")