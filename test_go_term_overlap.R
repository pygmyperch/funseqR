# Test GO Term Overlap Between Candidates and Background
# This script specifically tests the GO term overlap issue

cat("=== GO Term Overlap Analysis ===\n")

# Assume you have: annotations, candidates$bed_file

# Step 1: Get candidate locus IDs (replicate enrichment logic)
bed_coords <- candidates$bed_file
chrom_col <- if ("chrom" %in% colnames(bed_coords)) "chrom" else "chromosome"

candidate_df <- data.frame(
  chromosome = bed_coords[[chrom_col]],
  position = bed_coords$start + 1,
  stringsAsFactors = FALSE
)

candidate_patterns <- with(candidate_df, paste0("_", chromosome, "_", position, "$"))

candidate_locus_ids <- annotations$locus_id[
  sapply(annotations$locus_id, function(locus_id) {
    any(sapply(candidate_patterns, function(pattern) {
      grepl(pattern, locus_id)
    }))
  })
]

cat("Candidate locus IDs found:", length(candidate_locus_ids), "\n")

# Step 2: Get candidate annotations with GO terms
candidate_annotations <- annotations[annotations$locus_id %in% candidate_locus_ids, ]
candidate_with_go <- candidate_annotations[!is.na(candidate_annotations$go_terms) & 
                                          candidate_annotations$go_terms != "", ]

cat("Candidates with GO annotations:", nrow(candidate_with_go), "\n")

if (nrow(candidate_with_go) > 0) {
  cat("\nCandidate GO data (first 3):\n")
  for (i in 1:min(3, nrow(candidate_with_go))) {
    cat("Locus:", candidate_with_go$locus_id[i], "\n")
    cat("  GO terms:", substr(candidate_with_go$go_terms[i], 1, 80), "...\n")
    cat("  GO categories:", candidate_with_go$go_categories[i], "\n")
  }
  
  # Step 3: Extract all GO terms from candidates
  all_candidate_go_terms <- character(0)
  all_candidate_categories <- character(0)
  
  for (i in 1:nrow(candidate_with_go)) {
    terms <- unlist(strsplit(candidate_with_go$go_terms[i], ";"))
    categories <- unlist(strsplit(candidate_with_go$go_categories[i], ";"))
    
    # Clean terms
    terms <- terms[!is.na(terms) & terms != ""]
    categories <- categories[!is.na(categories) & categories != ""]
    
    all_candidate_go_terms <- c(all_candidate_go_terms, terms)
    all_candidate_categories <- c(all_candidate_categories, categories)
  }
  
  cat("\nCandidate GO summary:\n")
  cat("  Total GO term instances:", length(all_candidate_go_terms), "\n")
  cat("  Unique GO terms:", length(unique(all_candidate_go_terms)), "\n")
  
  # Count by category
  if (length(all_candidate_categories) > 0) {
    cat("  Categories:", paste(names(table(all_candidate_categories)), collapse = ", "), "\n")
    cat("  Category counts:", paste(table(all_candidate_categories), collapse = ", "), "\n")
  }
  
  # Step 4: Find most common candidate GO terms
  candidate_term_counts <- table(all_candidate_go_terms)
  top_candidate_terms <- head(sort(candidate_term_counts, decreasing = TRUE), 10)
  
  cat("\nTop 10 GO terms in candidates:\n")
  for (i in 1:length(top_candidate_terms)) {
    term_id <- names(top_candidate_terms)[i]
    candidate_count <- top_candidate_terms[i]
    
    # Count in full background
    background_matches <- grepl(paste0("(^|;)", term_id, "(;|$)"), annotations$go_terms)
    background_count <- sum(background_matches, na.rm = TRUE)
    
    cat(sprintf("  %s: %d candidates, %d background\n", term_id, candidate_count, background_count))
    
    # Calculate basic enrichment ratio
    if (background_count > 0) {
      candidate_freq <- candidate_count / length(candidate_locus_ids)
      background_freq <- background_count / nrow(annotations)
      enrichment_ratio <- candidate_freq / background_freq
      cat(sprintf("    Enrichment ratio: %.2f\n", enrichment_ratio))
    }
  }
  
  # Step 5: Focus on BP terms
  cat("\nBP (Biological Process) terms analysis:\n")
  
  # Get BP terms from candidates
  bp_indices <- which(all_candidate_categories == "P")
  if (length(bp_indices) > 0) {
    bp_terms_in_candidates <- all_candidate_go_terms[bp_indices]
    bp_term_counts <- table(bp_terms_in_candidates)
    
    cat("  BP terms in candidates:", length(unique(bp_terms_in_candidates)), "\n")
    cat("  Top BP terms:\n")
    
    top_bp <- head(sort(bp_term_counts, decreasing = TRUE), 5)
    for (i in 1:length(top_bp)) {
      term_id <- names(top_bp)[i]
      candidate_count <- top_bp[i]
      
      # Background count for this specific BP term
      background_matches <- grepl(paste0("(^|;)", term_id, "(;|$)"), annotations$go_terms)
      background_count <- sum(background_matches, na.rm = TRUE)
      
      cat(sprintf("    %s: %d candidates, %d background\n", term_id, candidate_count, background_count))
    }
  } else {
    cat("  No BP terms found in candidates\n")
  }
  
} else {
  cat("ERROR: No candidates have GO annotations!\n")
}

cat("\n=== Analysis Complete ===\n")