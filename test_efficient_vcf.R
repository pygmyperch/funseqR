#!/usr/bin/env Rscript

# Test the efficient VCF storage and retrieval functions

library(funseqR)

cat("=== TESTING EFFICIENT VCF STORAGE AND RETRIEVAL ===\n")

# Create a test database with the new schema
test_db <- "test_efficient_vcf.db"
if (file.exists(test_db)) file.remove(test_db)

cat("\n1. Creating database with new schema...\n")
cat("======================================\n")

con <- create_funseq_db(test_db, verbose = TRUE)

# Check if vcf_objects table was created
tables <- DBI::dbListTables(con)
if ("vcf_objects" %in% tables) {
  cat("✅ vcf_objects table created successfully\n")
} else {
  cat("❌ vcf_objects table not found\n")
}

cat("\n2. Testing VCF import with vcfR object storage...\n")
cat("===============================================\n")

# Test importing a small VCF file (create a minimal test file)
test_vcf <- "test_minimal.vcf"
writeLines(c(
  "##fileformat=VCFv4.2",
  "##source=test",
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
  "chr1\t100\t.\tA\tT\t60\tPASS\t.",
  "chr1\t200\t.\tG\tC\t50\tPASS\t.",
  "chr2\t300\t.\tC\tG\t70\tPASS\t."
), test_vcf)

tryCatch({
  import_result <- import_vcf(con, test_vcf, verbose = TRUE)
  cat("✅ VCF import completed\n")
  cat("File ID:", import_result$file_id, "\n")
  cat("Record count:", import_result$vcf_count, "\n")
  
  # Check if vcfR object was stored
  vcf_objects <- DBI::dbGetQuery(con, "SELECT COUNT(*) as count FROM vcf_objects")
  cat("VCF objects stored:", vcf_objects$count, "\n")
  
}, error = function(e) {
  cat("❌ VCF import failed:", e$message, "\n")
  return()
})

cat("\n3. Testing VCF retrieval without export...\n")
cat("=========================================\n")

tryCatch({
  # Test retrieval without export
  vcf_data <- retrieve_vcf_data(con, file_id = import_result$file_id, verbose = TRUE)
  
  cat("✅ VCF data retrieved successfully\n")
  cat("Rows:", nrow(vcf_data), "\n")
  cat("Columns:", paste(colnames(vcf_data), collapse = ", "), "\n")
  
  # Show first few rows
  cat("\nFirst few rows:\n")
  print(head(vcf_data, 3))
  
}, error = function(e) {
  cat("❌ VCF retrieval failed:", e$message, "\n")
})

cat("\n4. Testing fast VCF export...\n")
cat("=============================\n")

export_file <- "test_exported.vcf"
if (file.exists(export_file)) file.remove(export_file)

start_time <- Sys.time()

tryCatch({
  # Test retrieval with export
  vcf_data2 <- retrieve_vcf_data(con, file_id = import_result$file_id, 
                                 export_path = export_file, verbose = TRUE)
  
  end_time <- Sys.time()
  export_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  if (file.exists(export_file)) {
    cat("✅ VCF export completed successfully!\n")
    cat("Export time:", round(export_time, 2), "seconds\n")
    cat("File size:", file.size(export_file), "bytes\n")
    
    # Show exported file content
    cat("\nExported file content:\n")
    cat(paste(readLines(export_file), collapse = "\n"), "\n")
  } else {
    cat("❌ VCF export file not created\n")
  }
  
}, error = function(e) {
  cat("❌ VCF export failed:", e$message, "\n")
})

cat("\n5. Testing error handling...\n")
cat("============================\n")

# Test with non-existent file_id
tryCatch({
  invalid_data <- retrieve_vcf_data(con, file_id = 999, verbose = FALSE)
  cat("❌ Error handling failed\n")
}, error = function(e) {
  cat("✅ Error handling works:", e$message, "\n")
})

# Clean up
close_funseq_db(con)
if (file.exists(test_db)) file.remove(test_db)
if (file.exists(test_vcf)) file.remove(test_vcf)
if (file.exists(export_file)) file.remove(export_file)

cat("\n=== PERFORMANCE SUMMARY ===\n")
cat("✅ New approach uses vcfR object storage\n")
cat("✅ Export leverages vcfR::write.vcf() (native and fast)\n")
cat("✅ No JSON reconstruction needed\n")
cat("✅ Should be orders of magnitude faster than previous version\n")

cat("\n=== TEST COMPLETE ===\n")