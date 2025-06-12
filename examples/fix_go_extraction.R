# Quick Fix for Missing GO Terms Issue
#
# This script provides immediate solutions for the GO terms extraction problem

library(funseqR)

cat("=== GO Terms Extraction Issue - Quick Fixes ===\n\n")

# SOLUTION 1: Include IEA evidence codes
cat("SOLUTION 1: Include IEA Evidence Codes\n")
cat("========================================\n")
cat("The most likely cause is that your GO terms have IEA evidence codes,\n")
cat("which are excluded by default. IEA (Inferred from Electronic Annotation)\n")
cat("is very common in automated annotations.\n\n")

cat("Current default evidence_keep filter:\n")
cat("c('EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC')\n\n")

cat("Try re-running annotation with IEA included:\n\n")
cat("annotation_results <- annotate_blast_results(con,\n")
cat("                         blast_results$blast_param_id,\n")
cat("                         evidence_keep = c('EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC', 'IEA'),\n")
cat("                         verbose = TRUE)\n\n")

# SOLUTION 2: Debug specific accessions
cat("SOLUTION 2: Debug Specific Accessions\n")
cat("=====================================\n")
cat("Test a few specific accessions to see what's happening:\n\n")
cat("annotation_results <- annotate_blast_results(con,\n")
cat("                         blast_results$blast_param_id,\n")
cat("                         debug_accessions = c('Q49GP3', 'B2Z449', 'E7F4Z4'),\n")
cat("                         evidence_keep = c('EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC', 'IEA'),\n")
cat("                         verbose = TRUE)\n\n")

# SOLUTION 3: Check UniProt cache
cat("SOLUTION 3: Check UniProt Cache\n")
cat("===============================\n")
cat("Check what's stored in your UniProt cache:\n\n")
cat("# Check cache contents\n")
cat("cache_sample <- DBI::dbGetQuery(con, 'SELECT accession, data FROM uniprot_cache LIMIT 5')\n")
cat("print(cache_sample)\n\n")

cat("# Check for GO terms in cached data\n")
cat("for (i in 1:nrow(cache_sample)) {\n")
cat("  json_data <- jsonlite::fromJSON(cache_sample$data[i], flatten = FALSE)\n")
cat("  go_count <- length(grep('GO:', json_data, ignore.case = TRUE))\n")
cat("  cat('Accession', cache_sample$accession[i], 'has', go_count, 'GO references\\n')\n")
cat("}\n\n")

# SOLUTION 4: Completely disable evidence filtering
cat("SOLUTION 4: Disable Evidence Filtering\n")
cat("======================================\n")
cat("If you want ALL GO terms regardless of evidence:\n\n")
cat("annotation_results <- annotate_blast_results(con,\n")
cat("                         blast_results$blast_param_id,\n")
cat("                         evidence_keep = NULL,  # This disables evidence filtering\n")
cat("                         verbose = TRUE)\n\n")

# Quick diagnostic function
cat("DIAGNOSTIC FUNCTION\n")
cat("==================\n")
cat("Run this to quickly diagnose the issue:\n\n")

diagnose_go_issue <- function(con, blast_param_id) {
  cat("=== Diagnosing GO Terms Issue ===\n")
  
  # Get a sample of accessions
  accessions <- DBI::dbGetQuery(con, 
    "SELECT DISTINCT hit_accession FROM blast_results WHERE blast_param_id = ? LIMIT 5",
    params = list(blast_param_id))$hit_accession
  
  cat("Testing", length(accessions), "sample accessions:\n")
  
  for (acc in accessions) {
    cat("\nTesting", acc, ":\n")
    
    # Check cache
    cached <- DBI::dbGetQuery(con,
      "SELECT data FROM uniprot_cache WHERE accession = ?",
      params = list(acc))
    
    if (nrow(cached) > 0) {
      cat("  Found in cache\n")
      
      # Check for GO terms in cached data
      go_count <- length(grep("GO:", cached$data[1]))
      iea_count <- length(grep("IEA", cached$data[1]))
      
      cat("  GO references in raw data:", go_count, "\n")
      cat("  IEA evidence codes:", iea_count, "\n")
      
      # Test extraction
      tryCatch({
        info <- process_uniprot_json(cached$data[1], debug = FALSE)
        cat("  Extracted GO terms:", nrow(info$go_terms), "\n")
        if (nrow(info$go_terms) > 0) {
          evidence_types <- table(info$go_terms$go_evidence)
          cat("  Evidence types:", paste(names(evidence_types), evidence_types, collapse = ", "), "\n")
        }
      }, error = function(e) {
        cat("  Error in extraction:", e$message, "\n")
      })
    } else {
      cat("  Not in cache - would need API call\n")
    }
  }
  
  cat("\n=== Recommendations ===\n")
  if (any(grepl("IEA", accessions))) {
    cat("- Include IEA evidence codes in evidence_keep parameter\n")
  }
  cat("- Try verbose=TRUE to see detailed extraction process\n")
  cat("- Check if UniProt API responses have changed format\n")
}

cat("# Run the diagnostic\n")
cat("diagnose_go_issue(con, blast_results$blast_param_id)\n\n")

cat("=== Most Likely Solution ===\n")
cat("Based on the pattern (830 successful extractions but 0 GO terms),\n")
cat("the issue is almost certainly the evidence code filter.\n")
cat("Try this immediately:\n\n")
cat("annotation_results <- annotate_blast_results(con,\n")
cat("                         blast_results$blast_param_id,\n")
cat("                         evidence_keep = c('EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC', 'IEA'),\n")
cat("                         verbose = TRUE)\n\n")
cat("This should resolve the issue immediately!\n")