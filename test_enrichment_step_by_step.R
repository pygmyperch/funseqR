# Test GO Enrichment Step by Step
# This script tests each step of the GO enrichment process to identify issues

cat("=== GO Enrichment Step-by-Step Test ===\n")

# Assuming you have already run:
# annotations <- process_annotations(con, include = c("GO", "KEGG"))
# candidates <- import_candidate_loci(con, candidate_vcf_file, background_file_id = 1)

# TEST 1: Check annotation data structure
cat("\n1. ANNOTATION DATA STRUCTURE:\n")
cat("  - Total annotations:", nrow(annotations), "\n")
cat("  - Columns:", paste(colnames(annotations), collapse = ", "), "\n")

# Check GO annotation presence
go_present <- !is.na(annotations$go_terms) & annotations$go_terms != ""
cat("  - Annotations with GO terms:", sum(go_present), "\n")

if (sum(go_present) > 0) {
  # Show example GO data
  first_go <- annotations[go_present, ][1, ]
  cat("  - Example GO terms:", substr(first_go$go_terms, 1, 100), "...\n")
  cat("  - Example GO categories:", substr(first_go$go_categories, 1, 50), "...\n")
}

# TEST 2: Check candidate loci data structure  
cat("\n2. CANDIDATE LOCI STRUCTURE:\n")
cat("  - Type:", class(candidates$bed_file), "\n")
cat("  - Dimensions:", paste(dim(candidates$bed_file), collapse = " x "), "\n")
cat("  - Columns:", paste(colnames(candidates$bed_file), collapse = ", "), "\n")
cat("  - First few rows:\n")
print(head(candidates$bed_file, 3))

# TEST 3: Manual coordinate matching
cat("\n3. COORDINATE MATCHING TEST:\n")

# Extract coordinates from BED file
bed_coords <- candidates$bed_file
chrom_col <- if ("chrom" %in% colnames(bed_coords)) "chrom" else "chromosome"

# Convert BED to VCF-like coordinates (add 1 to start)
candidate_df <- data.frame(
  chromosome = bed_coords[[chrom_col]],
  position = bed_coords$start + 1,
  stringsAsFactors = FALSE
)

cat("  - Candidate coordinates (first 3):\n")
print(head(candidate_df, 3))

# Create matching patterns
candidate_patterns <- with(candidate_df, paste0("_", chromosome, "_", position, "$"))
cat("  - First few patterns:", paste(head(candidate_patterns, 3), collapse = ", "), "\n")

# Show some annotation locus_ids for comparison
cat("  - First few annotation locus_ids:", paste(head(annotations$locus_id, 3), collapse = ", "), "\n")

# Perform matching
matches <- sapply(annotations$locus_id, function(locus_id) {
  any(sapply(candidate_patterns, function(pattern) {
    grepl(pattern, locus_id)
  }))
})

candidate_locus_ids <- annotations$locus_id[matches]
cat("  - Matched locus IDs:", length(candidate_locus_ids), "\n")

if (length(candidate_locus_ids) > 0) {
  cat("  - First few matches:", paste(head(candidate_locus_ids, 3), collapse = ", "), "\n")
} else {
  cat("  - ERROR: No coordinate matches found!\n")
  cat("  - This explains why enrichment fails\n")
}

# TEST 4: GO category analysis
cat("\n4. GO CATEGORY ANALYSIS:\n")

if (sum(go_present) > 0) {
  # Extract all GO categories
  all_go_categories <- unlist(strsplit(annotations$go_categories[go_present], ";"))
  all_go_categories <- all_go_categories[!is.na(all_go_categories) & all_go_categories != ""]
  
  category_counts <- table(all_go_categories)
  cat("  - GO category distribution:\n")
  for (cat_name in names(category_counts)) {
    cat("    ", cat_name, ": ", category_counts[cat_name], "\n")
  }
  
  # Check if we have the expected single-letter codes
  expected_codes <- c("P", "F", "C")
  found_codes <- intersect(names(category_counts), expected_codes)
  missing_codes <- setdiff(expected_codes, names(category_counts))
  
  if (length(found_codes) > 0) {
    cat("  - Found expected codes:", paste(found_codes, collapse = ", "), "\n")
  }
  if (length(missing_codes) > 0) {
    cat("  - Missing expected codes:", paste(missing_codes, collapse = ", "), "\n")
  }
} else {
  cat("  - No GO annotations to analyze\n")
}

# TEST 5: Simulate enrichment for one term
cat("\n5. ENRICHMENT SIMULATION:\n")

if (length(candidate_locus_ids) > 0 && sum(go_present) > 0) {
  
  # Get GO data for candidates
  candidate_go_data <- annotations[annotations$locus_id %in% candidate_locus_ids & go_present, ]
  cat("  - Candidates with GO annotations:", nrow(candidate_go_data), "\n")
  
  if (nrow(candidate_go_data) > 0) {
    # Extract all GO terms from candidates
    candidate_go_terms <- unlist(strsplit(candidate_go_data$go_terms, ";"))
    candidate_go_terms <- candidate_go_terms[!is.na(candidate_go_terms) & candidate_go_terms != ""]
    
    # Count most common terms in candidates
    term_counts <- table(candidate_go_terms)
    top_terms <- head(sort(term_counts, decreasing = TRUE), 5)
    
    cat("  - Most common GO terms in candidates:\n")
    for (i in 1:length(top_terms)) {
      term_id <- names(top_terms)[i]
      count <- top_terms[i]
      cat("    ", term_id, ": ", count, " candidates\n")
      
      # Calculate background frequency
      background_count <- sum(grepl(term_id, annotations$go_terms, fixed = TRUE))
      cat("      (", background_count, " total in background)\n")
    }
  }
  
} else {
  if (length(candidate_locus_ids) == 0) {
    cat("  - Cannot simulate: No candidate loci matched\n")
  }
  if (sum(go_present) == 0) {
    cat("  - Cannot simulate: No GO annotations present\n")
  }
}

cat("\n=== Test Complete ===\n")
cat("Run this script to diagnose your enrichment issues:\n")
cat("source('test_enrichment_step_by_step.R')\n")