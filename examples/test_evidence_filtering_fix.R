# Test the Fixed Evidence Filtering Logic
#
# This script demonstrates the fix for compound evidence codes in GO term filtering

library(funseqR)

cat("=== Evidence Filtering Fix Test ===\n\n")

# The problem: Compound evidence codes
cat("PROBLEM: Compound Evidence Codes\n")
cat("================================\n")
cat("Your UniProt data contains compound evidence codes like:\n")
cat("- IEA:UniProtKB-SubCell\n")
cat("- IBA:GO_Central\n")
cat("- ISS:UniProtKB\n")
cat("- IMP:ZFIN\n\n")

cat("But the old filtering logic was looking for exact matches:\n")
cat("evidence_keep = c('IEA', 'IMP', ...)\n")
cat("'IEA:UniProtKB-SubCell' %in% c('IEA', 'IMP') = FALSE  # ❌ FAILED\n\n")

# The solution: Extract primary evidence codes
cat("SOLUTION: Extract Primary Evidence Codes\n")
cat("========================================\n")
cat("The fix extracts the primary code (before the colon):\n")
cat("'IEA:UniProtKB-SubCell' -> 'IEA' -> TRUE  # ✅ PASSES\n")
cat("'IMP:ZFIN' -> 'IMP' -> TRUE              # ✅ PASSES\n")
cat("'IBA:GO_Central' -> 'IBA' -> FALSE       # ❌ FILTERED (IBA not in keep list)\n\n")

# Test the filtering logic
test_evidence_codes <- c(
  "IEA:UniProtKB-SubCell",
  "IBA:GO_Central", 
  "ISS:UniProtKB",
  "IMP:ZFIN",
  "IDA:UniProtKB",
  "IPI:UniProtKB",
  "EXP",
  "TAS"
)

evidence_keep_conservative <- c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "TAS", "IC")
evidence_keep_with_iea <- c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "TAS", "IC", "IEA")

cat("TESTING WITH CONSERVATIVE FILTER (excludes IEA):\n")
cat("evidence_keep = c('EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC')\n\n")

for (code in test_evidence_codes) {
  primary_evidence <- strsplit(code, ":")[[1]][1]
  passes_filter <- primary_evidence %in% evidence_keep_conservative
  status <- if (passes_filter) "✅ PASS" else "❌ FILTER"
  cat(sprintf("%-25s -> %-5s -> %s\n", code, primary_evidence, status))
}

cat("\nTESTING WITH IEA INCLUDED:\n")
cat("evidence_keep = c('EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC', 'IEA')\n\n")

for (code in test_evidence_codes) {
  primary_evidence <- strsplit(code, ":")[[1]][1]
  passes_filter <- primary_evidence %in% evidence_keep_with_iea
  status <- if (passes_filter) "✅ PASS" else "❌ FILTER"
  cat(sprintf("%-25s -> %-5s -> %s\n", code, primary_evidence, status))
}

cat("\n=== Now Test with Your Data ===\n")
cat("With the fix applied, you can now use:\n\n")

cat("# Conservative filtering (high-confidence only)\n")
cat("annotation_results <- annotate_blast_results(con,\n")
cat("                         blast_results$blast_param_id,\n")
cat("                         evidence_keep = c('EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC'),\n")
cat("                         verbose = TRUE)\n\n")

cat("# Include computational annotations (recommended for non-model organisms)\n")
cat("annotation_results <- annotate_blast_results(con,\n")
cat("                         blast_results$blast_param_id,\n")
cat("                         evidence_keep = c('EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC', 'IEA'),\n")
cat("                         verbose = TRUE)\n\n")

cat("# Accept all evidence types\n")
cat("annotation_results <- annotate_blast_results(con,\n")
cat("                         blast_results$blast_param_id,\n")
cat("                         evidence_keep = NULL,\n")
cat("                         verbose = TRUE)\n\n")

cat("=== Expected Results ===\n")
cat("With the fix, you should now get GO terms even with conservative filtering!\n")
cat("The compound evidence codes will be properly parsed and filtered.\n")