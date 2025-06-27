#!/usr/bin/env Rscript

# Investigation script for UniProt cache data
tryCatch({
  library(DBI)
  library(RSQLite) 
  library(jsonlite)
}, error = function(e) {
  cat("Error loading libraries:", e$message, "\n")
  quit(status = 1)
})

# Connect to database
tryCatch({
  con <- dbConnect(SQLite(), "funseq_project.db")
}, error = function(e) {
  cat("Error connecting to database:", e$message, "\n")
  quit(status = 1)
})

cat("=== UniProt Cache Investigation ===\n\n")

# Check available tables
tables <- dbListTables(con)
cat("Available tables:", paste(tables, collapse=", "), "\n\n")

if ("uniprot_cache" %in% tables) {
  # Get basic stats
  count <- dbGetQuery(con, "SELECT COUNT(*) as count FROM uniprot_cache")$count
  cat("Total UniProt cache entries:", count, "\n\n")
  
  if (count > 0) {
    # Get sample accessions and metadata
    sample_meta <- dbGetQuery(con, "
      SELECT accession, LENGTH(response_json) as json_length, retrieval_date 
      FROM uniprot_cache 
      ORDER BY json_length DESC 
      LIMIT 5
    ")
    
    cat("Sample cache entries (ordered by JSON size):\n")
    print(sample_meta)
    cat("\n")
    
    # Get one detailed response for analysis
    sample_response <- dbGetQuery(con, "
      SELECT accession, response_json 
      FROM uniprot_cache 
      WHERE LENGTH(response_json) > 1000 
      LIMIT 1
    ")
    
    if (nrow(sample_response) > 0) {
      accession <- sample_response$accession[1]
      json_text <- sample_response$response_json[1]
      
      cat("=== Detailed Analysis of", accession, "===\n")
      cat("JSON length:", nchar(json_text), "characters\n\n")
      
      # Parse JSON and analyze structure
      tryCatch({
        parsed <- fromJSON(json_text, simplifyVector = FALSE)
        
        # Find the entry
        entry <- NULL
        if (!is.null(parsed$results) && length(parsed$results) > 0) {
          entry <- parsed$results[[1]]
        } else if (!is.null(parsed$primaryAccession)) {
          entry <- parsed
        }
        
        if (!is.null(entry)) {
          cat("Entry found for accession:", entry$primaryAccession %||% "Unknown", "\n")
          cat("UniProt ID:", entry$uniProtkbId %||% "Unknown", "\n\n")
          
          # Analyze cross-references
          if (!is.null(entry$uniProtKBCrossReferences)) {
            cat("=== Cross-References Analysis ===\n")
            
            # Count cross-references by database
            db_counts <- table(sapply(entry$uniProtKBCrossReferences, function(x) x$database %||% "Unknown"))
            cat("Cross-reference databases found:\n")
            for (db in names(db_counts)) {
              cat("  -", db, ":", db_counts[db], "entries\n")
            }
            cat("\n")
            
            # Examine specific annotation types we're interested in
            target_dbs <- c("Pfam", "InterPro", "eggNOG", "SubCell", "HAMAP", "SUPFAM", "PROSITE", "CDD")
            
            for (target_db in target_dbs) {
              refs <- entry$uniProtKBCrossReferences[sapply(entry$uniProtKBCrossReferences, function(x) x$database == target_db)]
              
              if (length(refs) > 0) {
                cat("=== ", target_db, " References ===\n")
                for (i in seq_along(refs)) {
                  ref <- refs[[i]]
                  cat("ID:", ref$id %||% "Unknown", "\n")
                  
                  if (!is.null(ref$properties) && length(ref$properties) > 0) {
                    cat("Properties:\n")
                    for (prop in ref$properties) {
                      cat("  ", prop$key, ":", prop$value, "\n")
                    }
                  }
                  cat("\n")
                }
              }
            }
            
            # Show first few cross-references for detailed structure
            cat("=== First 3 Cross-References (detailed structure) ===\n")
            for (i in 1:min(3, length(entry$uniProtKBCrossReferences))) {
              ref <- entry$uniProtKBCrossReferences[[i]]
              cat("Reference", i, ":\n")
              cat("  Database:", ref$database %||% "Unknown", "\n")
              cat("  ID:", ref$id %||% "Unknown", "\n")
              if (!is.null(ref$properties)) {
                cat("  Properties:\n")
                for (prop in ref$properties) {
                  cat("    ", prop$key, ":", prop$value, "\n")
                }
              }
              cat("\n")
            }
            
          } else {
            cat("No cross-references found in this entry\n")
          }
          
        } else {
          cat("Could not find entry data in JSON\n")
        }
        
      }, error = function(e) {
        cat("Error parsing JSON:", e$message, "\n")
        cat("First 500 characters of JSON:\n")
        cat(substr(json_text, 1, 500), "\n")
      })
    }
  }
} else {
  cat("uniprot_cache table not found\n")
}

dbDisconnect(con)
cat("\nInvestigation complete.\n")