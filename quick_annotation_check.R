#!/usr/bin/env Rscript
#
# Quick UniProt Annotation Check
# Simple script to verify annotation availability
#

library(DBI)
library(RSQLite)
library(jsonlite)

# Quick check function
quick_check <- function() {
  cat("=== Quick UniProt Annotation Check ===\n\n")
  
  con <- dbConnect(SQLite(), "funseq_project.db")
  on.exit(dbDisconnect(con))
  
  # Get one good example
  sample <- dbGetQuery(con, "
    SELECT accession, response_json 
    FROM uniprot_cache 
    WHERE LENGTH(response_json) > 2000
    LIMIT 1
  ")
  
  if (nrow(sample) == 0) {
    cat("No suitable cache entries found\n")
    return()
  }
  
  cat("Analyzing:", sample$accession, "\n\n")
  
  # Parse and examine structure
  parsed <- fromJSON(sample$response_json, simplifyVector = FALSE)
  
  # Find entry
  entry <- NULL
  if (!is.null(parsed$results) && length(parsed$results) > 0) {
    entry <- parsed$results[[1]]
  } else if (!is.null(parsed$primaryAccession)) {
    entry <- parsed
  }
  
  if (is.null(entry$uniProtKBCrossReferences)) {
    cat("No cross-references found\n")
    return()
  }
  
  # Count by database
  databases <- table(sapply(entry$uniProtKBCrossReferences, function(x) x$database))
  
  cat("Cross-reference databases found:\n")
  for (db in names(databases)) {
    cat("  ", db, ":", databases[db], "\n")
  }
  
  cat("\n=== Key Findings ===\n")
  
  # Check for our target annotation types
  targets <- c("Pfam", "InterPro", "eggNOG", "SubCell")
  found <- intersect(names(databases), targets)
  missing <- setdiff(targets, names(databases))
  
  if (length(found) > 0) {
    cat("✅ Available annotation types:\n")
    for (f in found) {
      cat("   ", f, "\n")
    }
  }
  
  if (length(missing) > 0) {
    cat("❌ Missing annotation types:\n")
    for (m in missing) {
      cat("   ", m, "\n")
    }
  }
  
  # Show example structure for first target found
  if (length(found) > 0) {
    target_db <- found[1]
    cat("\n=== Example", target_db, "structure ===\n")
    
    target_refs <- entry$uniProtKBCrossReferences[sapply(entry$uniProtKBCrossReferences, function(x) x$database == target_db)]
    if (length(target_refs) > 0) {
      ref <- target_refs[[1]]
      cat("ID:", ref$id, "\n")
      if (!is.null(ref$properties)) {
        cat("Properties:\n")
        for (prop in ref$properties) {
          cat("  ", prop$key, ":", prop$value, "\n")
        }
      }
    }
  }
}

# Run if called as script
if (!interactive()) {
  quick_check()
}