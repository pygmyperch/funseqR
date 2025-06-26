# Debug GO Enrichment Analysis
# This script helps diagnose why GO enrichment is not finding significant terms

# Function to debug GO enrichment step by step
debug_go_enrichment <- function(annotations, candidate_loci, ontology = "BP") {
  
  message("=== GO Enrichment Debug Analysis ===")
  message("Ontology: ", ontology)
  
  # Step 1: Check candidate loci identification
  message("\n1. Candidate Loci Identification:")
  
  # Simulate the coordinate matching logic from the enrichment function
  has_chrom <- "chrom" %in% colnames(candidate_loci)
  has_chromosome <- "chromosome" %in% colnames(candidate_loci)
  has_start <- "start" %in% colnames(candidate_loci)
  has_end <- "end" %in% colnames(candidate_loci)
  
  if ((has_chrom || has_chromosome) && has_start && has_end) {
    chrom_col <- if (has_chrom) "chrom" else "chromosome"
    message("  - Detected BED format: ", chrom_col, ", start, end")
    
    candidate_df <- data.frame(
      chromosome = candidate_loci[[chrom_col]],
      position = candidate_loci$start + 1,  # Convert 0-based to 1-based
      stringsAsFactors = FALSE
    )
  }
  
  # Create matching patterns
  candidate_patterns <- with(candidate_df, paste0("_", chromosome, "_", position, "$"))
  
  # Find matching loci
  candidate_locus_ids <- annotations$locus_id[
    sapply(annotations$locus_id, function(locus_id) {
      any(sapply(candidate_patterns, function(pattern) {
        grepl(pattern, locus_id)
      }))
    })
  ]
  
  message("  - Total candidate coordinates: ", nrow(candidate_df))
  message("  - First few patterns: ", paste(head(candidate_patterns, 3), collapse = ", "))
  message("  - Matched locus IDs: ", length(candidate_locus_ids))
  message("  - First few matches: ", paste(head(candidate_locus_ids, 3), collapse = ", "))
  
  # Step 2: Check GO data preparation
  message("\n2. GO Data Preparation:")
  
  # Filter annotations with GO terms
  go_annotations <- annotations[!is.na(annotations$go_terms) & annotations$go_terms != "", ]
  message("  - Total annotations: ", nrow(annotations))
  message("  - Annotations with GO terms: ", nrow(go_annotations))
  
  # Create GO term mappings (simplified version)
  go_term_map <- list()
  
  for (i in 1:min(nrow(go_annotations), 100)) {  # Limit to first 100 for debugging
    locus_id <- go_annotations$locus_id[i]
    go_terms <- unlist(strsplit(go_annotations$go_terms[i], ";"))
    go_categories <- unlist(strsplit(go_annotations$go_categories[i], ";"))
    
    for (j in seq_along(go_terms)) {
      go_id <- go_terms[j]
      if (!is.na(go_id) && go_id != "") {
        if (!go_id %in% names(go_term_map)) {
          go_term_map[[go_id]] <- list(
            go_id = go_id,
            go_category = if (j <= length(go_categories)) go_categories[j] else "Unknown",
            loci = character(0)
          )
        }
        go_term_map[[go_id]]$loci <- c(go_term_map[[go_id]]$loci, locus_id)
      }
    }
  }
  
  message("  - GO terms found: ", length(go_term_map))
  
  # Step 3: Check category mapping
  message("\n3. Category Analysis:")
  
  ontology_map <- c("BP" = "P", "MF" = "F", "CC" = "C")
  category_code <- ontology_map[ontology]
  
  message("  - Looking for ontology: ", ontology, " (code: ", category_code, ")")
  
  # Count terms by category
  all_categories <- sapply(go_term_map, function(x) x$go_category)
  category_counts <- table(all_categories)
  message("  - Available categories: ", paste(names(category_counts), "=", category_counts, collapse = ", "))
  
  # Filter for specific ontology
  ontology_terms <- go_term_map[sapply(go_term_map, function(x) x$go_category == category_code)]
  message("  - Terms in ", ontology, " ontology: ", length(ontology_terms))
  
  # Step 4: Check candidate/background overlap
  message("\n4. Candidate/Background Overlap:")
  
  candidate_with_go <- intersect(candidate_locus_ids, go_annotations$locus_id)
  message("  - Candidate loci with GO annotations: ", length(candidate_with_go))
  
  if (length(ontology_terms) > 0) {
    # Check first few terms for overlap
    for (i in 1:min(3, length(ontology_terms))) {
      term_info <- ontology_terms[[i]]
      term_loci <- term_info$loci
      candidate_overlap <- length(intersect(candidate_locus_ids, term_loci))
      background_total <- length(intersect(go_annotations$locus_id, term_loci))
      
      message("  - Term ", term_info$go_id, ": ", 
              candidate_overlap, " candidates / ", background_total, " background")
    }
  }
  
  # Step 5: Statistical testing preview
  message("\n5. Statistical Testing Feasibility:")
  
  if (length(candidate_locus_ids) > 0 && length(ontology_terms) > 0) {
    message("  - Candidate set size: ", length(candidate_locus_ids))
    message("  - Background set size: ", nrow(go_annotations))
    message("  - Terms to test: ", length(ontology_terms))
    
    # Check if we have enough candidates for meaningful testing
    if (length(candidate_locus_ids) < 5) {
      message("  - WARNING: Very small candidate set may limit statistical power")
    }
    
    if (length(candidate_with_go) == 0) {
      message("  - ERROR: No candidates have GO annotations - enrichment impossible")
    }
  }
  
  message("\n=== Debug Complete ===")
  
  return(list(
    candidate_locus_ids = candidate_locus_ids,
    candidate_with_go = candidate_with_go,
    go_term_count = length(go_term_map),
    ontology_term_count = length(ontology_terms),
    category_counts = category_counts
  ))
}

# Usage example:
# debug_results <- debug_go_enrichment(annotations, candidates$bed_file, "BP")