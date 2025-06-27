#!/usr/bin/env Rscript
#
# UniProt Cache Annotation Verification Script
# This script examines cached UniProt API responses to identify available annotation types
#

library(DBI)
library(RSQLite)
library(jsonlite)

# Function to safely extract nested values
safe_extract <- function(obj, path) {
  tryCatch({
    if (is.null(obj)) return(NA)
    result <- obj
    for (p in path) {
      if (is.null(result[[p]])) return(NA)
      result <- result[[p]]
    }
    return(result)
  }, error = function(e) NA)
}

# Main analysis function
analyze_uniprot_cache <- function(db_path = "funseq_project.db") {
  cat("=== UniProt Cache Annotation Analysis ===\n\n")
  
  # Connect to database
  con <- dbConnect(SQLite(), db_path)
  on.exit(dbDisconnect(con))
  
  # Check if cache table exists
  tables <- dbListTables(con)
  if (!"uniprot_cache" %in% tables) {
    cat("ERROR: uniprot_cache table not found in database\n")
    cat("Available tables:", paste(tables, collapse = ", "), "\n")
    return(invisible(NULL))
  }
  
  # Get cache statistics
  cache_stats <- dbGetQuery(con, "
    SELECT 
      COUNT(*) as total_entries,
      MIN(LENGTH(response_json)) as min_json_size,
      MAX(LENGTH(response_json)) as max_json_size,
      AVG(LENGTH(response_json)) as avg_json_size
    FROM uniprot_cache 
    WHERE response_json IS NOT NULL 
    AND response_json != ''
    AND response_json NOT LIKE '%\"empty\":true%'
  ")
  
  cat("Cache Statistics:\n")
  cat("  Total entries:", cache_stats$total_entries, "\n")
  cat("  JSON size range:", cache_stats$min_json_size, "-", cache_stats$max_json_size, "chars\n")
  cat("  Average JSON size:", round(cache_stats$avg_json_size), "chars\n\n")
  
  if (cache_stats$total_entries == 0) {
    cat("No valid cache entries found for analysis\n")
    return(invisible(NULL))
  }
  
  # Get sample entries for analysis
  sample_entries <- dbGetQuery(con, "
    SELECT accession, response_json 
    FROM uniprot_cache 
    WHERE response_json IS NOT NULL 
    AND response_json != ''
    AND response_json NOT LIKE '%\"empty\":true%'
    AND LENGTH(response_json) > 1000
    ORDER BY LENGTH(response_json) DESC
    LIMIT 10
  ")
  
  cat("Analyzing", nrow(sample_entries), "sample entries...\n\n")
  
  # Initialize tracking variables
  all_databases <- list()
  annotation_examples <- list()
  property_keys <- list()
  
  # Analyze each sample entry
  for (i in 1:nrow(sample_entries)) {
    accession <- sample_entries$accession[i]
    json_text <- sample_entries$response_json[i]
    
    cat("--- Analyzing", accession, "---\n")
    
    tryCatch({
      # Parse JSON
      parsed <- fromJSON(json_text, simplifyVector = FALSE)
      
      # Find entry data
      entry <- NULL
      if (!is.null(parsed$results) && length(parsed$results) > 0) {
        entry <- parsed$results[[1]]
      } else if (!is.null(parsed$primaryAccession)) {
        entry <- parsed
      }
      
      if (is.null(entry)) {
        cat("  Could not find entry data\n")
        next
      }
      
      # Analyze cross-references
      xrefs <- safe_extract(entry, "uniProtKBCrossReferences")
      if (is.null(xrefs) || length(xrefs) == 0) {
        cat("  No cross-references found\n")
        next
      }
      
      cat("  Found", length(xrefs), "cross-references\n")
      
      # Count databases
      entry_databases <- list()
      for (ref in xrefs) {
        db_name <- safe_extract(ref, "database")
        if (!is.na(db_name)) {
          if (is.null(entry_databases[[db_name]])) {
            entry_databases[[db_name]] <- 0
          }
          entry_databases[[db_name]] <- entry_databases[[db_name]] + 1
          
          # Store example for each database
          if (is.null(annotation_examples[[db_name]])) {
            annotation_examples[[db_name]] <- list(
              accession = accession,
              id = safe_extract(ref, "id"),
              properties = safe_extract(ref, "properties")
            )
          }
          
          # Collect property keys
          props <- safe_extract(ref, "properties")
          if (!is.null(props) && length(props) > 0) {
            for (prop in props) {
              key <- safe_extract(prop, "key")
              if (!is.na(key)) {
                if (is.null(property_keys[[db_name]])) {
                  property_keys[[db_name]] <- list()
                }
                property_keys[[db_name]][[key]] <- TRUE
              }
            }
          }
        }
      }
      
      # Update global database counts
      for (db in names(entry_databases)) {
        if (is.null(all_databases[[db]])) {
          all_databases[[db]] <- list(count = 0, entries = 0)
        }
        all_databases[[db]]$count <- all_databases[[db]]$count + entry_databases[[db]]
        all_databases[[db]]$entries <- all_databases[[db]]$entries + 1
      }
      
      # Show databases for this entry
      if (length(entry_databases) > 0) {
        db_summary <- paste(names(entry_databases), collapse = ", ")
        cat("  Databases:", db_summary, "\n")
      }
      
    }, error = function(e) {
      cat("  Error parsing JSON:", e$message, "\n")
    })
  }
  
  # Summary Report
  cat("\n=== ANNOTATION TYPES SUMMARY ===\n\n")
  
  if (length(all_databases) > 0) {
    # Sort databases by frequency
    db_freq <- sapply(all_databases, function(x) x$count)
    sorted_dbs <- names(sort(db_freq, decreasing = TRUE))
    
    cat("Available annotation databases (sorted by frequency):\n")
    for (db in sorted_dbs) {
      info <- all_databases[[db]]
      cat(sprintf("  %-12s: %3d annotations across %d entries\n", 
                  db, info$count, info$entries))
    }
    
    cat("\n=== DETAILED ANNOTATION STRUCTURE ===\n\n")
    
    # Focus on annotation types we're interested in
    target_dbs <- c("Pfam", "InterPro", "eggNOG", "SubCell", "HAMAP", "SUPFAM", "PROSITE", "CDD", "GO", "KEGG")
    available_targets <- intersect(target_dbs, names(all_databases))
    
    for (db in available_targets) {
      cat("=== ", db, " ===\n")
      
      example <- annotation_examples[[db]]
      if (!is.null(example)) {
        cat("Example from", example$accession, ":\n")
        cat("  ID:", example$id %||% "Unknown", "\n")
        
        if (!is.null(example$properties) && length(example$properties) > 0) {
          cat("  Properties:\n")
          for (prop in example$properties) {
            key <- safe_extract(prop, "key")
            value <- safe_extract(prop, "value")
            if (!is.na(key) && !is.na(value)) {
              cat("    ", key, ":", value, "\n")
            }
          }
        }
      }
      
      # Show all property keys found for this database
      if (!is.null(property_keys[[db]])) {
        all_keys <- names(property_keys[[db]])
        cat("  All property keys found:", paste(all_keys, collapse = ", "), "\n")
      }
      
      cat("\n")
    }
    
    # Recommendations
    cat("=== IMPLEMENTATION RECOMMENDATIONS ===\n\n")
    
    recommended <- intersect(c("Pfam", "InterPro", "eggNOG", "SubCell"), available_targets)
    if (length(recommended) > 0) {
      cat("✅ RECOMMENDED for implementation:\n")
      for (db in recommended) {
        info <- all_databases[[db]]
        cat("  ", db, "- Available in", info$entries, "entries\n")
      }
    }
    
    additional <- setdiff(available_targets, c("GO", "KEGG", recommended))
    if (length(additional) > 0) {
      cat("\n🔧 ADDITIONAL options:\n")
      for (db in additional) {
        info <- all_databases[[db]]
        cat("  ", db, "- Available in", info$entries, "entries\n")
      }
    }
    
    missing <- setdiff(target_dbs, available_targets)
    if (length(missing) > 0) {
      cat("\n❌ NOT FOUND in cache:\n")
      for (db in missing) {
        cat("  ", db, "\n")
      }
    }
    
  } else {
    cat("No cross-reference databases found in any entries\n")
  }
  
  cat("\n=== NEXT STEPS ===\n")
  cat("1. Review the available annotation types above\n")
  cat("2. Decide which additional types to implement\n")
  cat("3. Use the property key information to extend extract_uniprot_info()\n")
  cat("4. Add corresponding database tables for storage\n\n")
  
  invisible(list(
    databases = all_databases,
    examples = annotation_examples,
    properties = property_keys
  ))
}

# Run analysis if called as script
if (!interactive()) {
  analyze_uniprot_cache()
}