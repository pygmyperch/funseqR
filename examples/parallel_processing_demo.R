# Enhanced Parallel Processing Demo for funseqR
# 
# This example demonstrates the new parallel processing capabilities with
# enhanced dual progress bars in import_flanking_seqs_to_db() for faster ORF extraction

library(funseqR)

# Connect to database
con <- connect_funseq_db("demo_analysis.db")

# Example 1: Sequential processing (threads = 1)
cat("=== Sequential Processing ===\n")
system.time({
  result_sequential <- import_flanking_seqs_to_db(
    con = con,
    vcf_file_id = 1,
    genome_id = 1,
    flank_size = 500,
    translate_flanks = TRUE,
    threads = 1,
    verbose = TRUE
  )
})

cat("\nSequential results:\n")
print(result_sequential)

# Example 2: Parallel processing (threads = 4)
cat("\n=== Parallel Processing ===\n")
system.time({
  result_parallel <- import_flanking_seqs_to_db(
    con = con,
    vcf_file_id = 2,  # Different VCF file
    genome_id = 1,
    flank_size = 500,
    translate_flanks = TRUE,
    threads = 4,
    batch_size = 1000,
    verbose = TRUE
  )
})

cat("\nParallel results:\n")
print(result_parallel)

# Example 3: Parallel processing with custom parameters
cat("\n=== Parallel Processing with Custom Parameters ===\n")
system.time({
  result_custom <- import_flanking_seqs_to_db(
    con = con,
    vcf_file_id = 3,  # Another VCF file
    genome_id = 1,
    flank_size = 1000,    # Larger flanking regions
    translate_flanks = TRUE,
    orf_min_aa = 50,      # Longer minimum ORF length
    orf_return = "aa",    # Return amino acid sequences
    threads = 8,          # More threads
    batch_size = 500,     # Smaller batch size
    verbose = TRUE
  )
})

cat("\nCustom parallel results:\n")
print(result_custom)

# Performance comparison note
cat("\n=== Enhanced Progress Bar Features ===\n")
cat("- Dual progress bars for parallel processing:\n")
cat("  * Phase 1: Parallel ORF computation with thread activity\n")
cat("  * Phase 2: Sequential database writing with throughput metrics\n")
cat("- Single enhanced progress bar for sequential processing\n")
cat("- Real-time performance metrics (sequences/sec, records/sec)\n")
cat("- Processing time breakdown and efficiency analysis\n")
cat("- No more 'hanging' appearance during parallel ORF extraction\n")
cat("\n=== Performance Notes ===\n")
cat("- Parallel processing is automatically enabled when threads > 1 and translate_flanks = TRUE\n")
cat("- Performance benefits are most noticeable with large datasets (>1000 sequences)\n")
cat("- The optimal number of threads depends on your system's CPU cores\n")
cat("- Use threads = parallel::detectCores() - 1 for maximum performance\n")
cat("- Sequential database writes ensure data integrity\n")

# Clean up
close_funseq_db(con)