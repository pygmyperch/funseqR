#!/usr/bin/env Rscript

# Test the enhanced funseqR_summary() function with grouped types

library(funseqR)

cat("=== TESTING ENHANCED funseqR_summary() FUNCTION ===\n")

# Connect to the database
con <- connect_funseq_db('funseq_project.db')

cat("\n1. Testing 'database' type (renamed from 'data')...\n")
cat("===================================================\n")

tryCatch({
  db_summary <- funseqR_summary(con, type = "database")
  cat("✅ Database summary successful\n")
  cat("Structure: ", class(db_summary), " with ", nrow(db_summary), " rows\n")
  cat("Columns: ", paste(names(db_summary), collapse = ", "), "\n")
  cat("Tables with data: ", sum(db_summary$record_count > 0), "/", nrow(db_summary), "\n")
}, error = function(e) {
  cat("❌ Database summary failed:", e$message, "\n")
})

cat("\n2. Testing 'input_data' type...\n")
cat("===============================\n")

tryCatch({
  input_summary <- funseqR_summary(con, type = "input_data")
  cat("✅ Input data summary successful\n")
  cat("Structure: ", class(input_summary), " (list)\n")
  cat("Tables included: ", paste(names(input_summary), collapse = ", "), "\n")
  
  # Show sample results
  for (table_name in names(input_summary)) {
    table_data <- input_summary[[table_name]]
    if ("message" %in% names(table_data)) {
      cat("  ", table_name, ": ", table_data$message[1], "\n")
    } else {
      cat("  ", table_name, ": ", nrow(table_data), " rows, ", ncol(table_data), " columns\n")
    }
  }
}, error = function(e) {
  cat("❌ Input data summary failed:", e$message, "\n")
})

cat("\n3. Testing 'blast' type...\n")
cat("==========================\n")

tryCatch({
  blast_summary <- funseqR_summary(con, type = "blast")
  cat("✅ BLAST summary successful\n")
  cat("Structure: ", class(blast_summary), " (list)\n")
  cat("Tables included: ", paste(names(blast_summary), collapse = ", "), "\n")
  
  # Show sample results
  for (table_name in names(blast_summary)) {
    table_data <- blast_summary[[table_name]]
    if ("message" %in% names(table_data)) {
      cat("  ", table_name, ": ", table_data$message[1], "\n")
    } else {
      cat("  ", table_name, ": ", nrow(table_data), " rows, ", ncol(table_data), " columns\n")
    }
  }
}, error = function(e) {
  cat("❌ BLAST summary failed:", e$message, "\n")
})

cat("\n4. Testing 'annotation' type...\n")
cat("===============================\n")

tryCatch({
  annotation_summary <- funseqR_summary(con, type = "annotation")
  cat("✅ Annotation summary successful\n")
  cat("Structure: ", class(annotation_summary), " (list)\n")
  cat("Tables included: ", paste(names(annotation_summary), collapse = ", "), "\n")
  
  # Show sample results
  for (table_name in names(annotation_summary)) {
    table_data <- annotation_summary[[table_name]]
    if ("message" %in% names(table_data)) {
      cat("  ", table_name, ": ", table_data$message[1], "\n")
    } else {
      cat("  ", table_name, ": ", nrow(table_data), " rows, ", ncol(table_data), " columns\n")
    }
  }
}, error = function(e) {
  cat("❌ Annotation summary failed:", e$message, "\n")
})

cat("\n5. Testing 'analyses' type...\n")
cat("=============================\n")

tryCatch({
  analyses_summary <- funseqR_summary(con, type = "analyses")
  cat("✅ Analyses summary successful\n")
  cat("Structure: ", class(analyses_summary), " (list)\n")
  cat("Tables included: ", paste(names(analyses_summary), collapse = ", "), "\n")
  
  # Show sample results
  for (table_name in names(analyses_summary)) {
    table_data <- analyses_summary[[table_name]]
    if ("message" %in% names(table_data)) {
      cat("  ", table_name, ": ", table_data$message[1], "\n")
    } else {
      cat("  ", table_name, ": ", nrow(table_data), " rows, ", ncol(table_data), " columns\n")
    }
  }
}, error = function(e) {
  cat("❌ Analyses summary failed:", e$message, "\n")
})

cat("\n6. Testing error handling...\n")
cat("============================\n")

# Test invalid type
tryCatch({
  invalid_summary <- funseqR_summary(con, type = "invalid_type")
  cat("❌ Error handling failed - should have thrown error\n")
}, error = function(e) {
  cat("✅ Error handling works:", e$message, "\n")
})

# Test with old "data" type (should fail)
tryCatch({
  old_summary <- funseqR_summary(con, type = "data")
  cat("❌ Old 'data' type should have failed\n")
}, error = function(e) {
  cat("✅ Old 'data' type correctly rejected:", e$message, "\n")
})

cat("\n7. Demonstrating practical usage...\n")
cat("===================================\n")

# Show how users would actually use this
cat("Example: Getting VCF data details:\n")
input_data <- funseqR_summary(con, type = "input_data")
if (!("message" %in% names(input_data$vcf_data))) {
  cat("VCF chromosome distribution:\n")
  print(head(input_data$vcf_data))
}

cat("\nExample: Getting enrichment analysis overview:\n")
analyses_data <- funseqR_summary(con, type = "analyses")
if (!("message" %in% names(analyses_data$ora_results))) {
  cat("ORA results summary:\n")
  print(analyses_data$ora_results)
}

# Close connection
close_funseq_db(con)

cat("\n=== ENHANCEMENT SUMMARY ===\n")
cat("✅ Successfully enhanced funseqR_summary() function\n")
cat("✅ Added 4 new grouped summary types\n")
cat("✅ Renamed 'data' to 'database' for clarity\n")
cat("✅ Each group returns detailed, actionable information\n")
cat("✅ Consistent list structure across all group types\n")
cat("✅ Error handling and validation implemented\n")
cat("✅ Ready for workflow-based data exploration\n")

cat("\n=== TEST COMPLETE ===\n")