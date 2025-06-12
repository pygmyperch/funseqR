# Debug GO Terms Extraction Issue
# 
# This script helps diagnose why GO terms are not being extracted from UniProt responses

library(funseqR)

# Function to debug a specific UniProt accession
debug_uniprot_response <- function(accession, verbose = TRUE) {
  cat("=== Debugging UniProt Response for", accession, "===\n")
  
  # Test the actual API call
  url <- paste0("https://rest.uniprot.org/uniprotkb/search?query=accession:", 
                accession, "&format=json")
  
  if (verbose) cat("API URL:", url, "\n")
  
  # Make the API call
  tryCatch({
    response <- httr::GET(url, httr::add_headers(
      "User-Agent" = "funseqR/R API client",
      "Accept" = "application/json"
    ))
    
    if (httr::status_code(response) == 200) {
      content <- httr::content(response, as = "text", encoding = "UTF-8")
      json_data <- jsonlite::fromJSON(content, flatten = FALSE)
      
      if (verbose) {
        cat("API call successful, status:", httr::status_code(response), "\n")
        cat("Response contains", length(json_data$results), "results\n")
      }
      
      if (length(json_data$results) > 0) {
        entry <- json_data$results[[1]]
        
        # Check for cross-references
        if (!is.null(entry$uniProtKBCrossReferences)) {
          cat("Found", length(entry$uniProtKBCrossReferences), "cross-references\n")
          
          # Count GO terms
          go_refs <- 0
          kegg_refs <- 0
          other_refs <- 0
          
          for (i in seq_along(entry$uniProtKBCrossReferences)) {
            ref <- entry$uniProtKBCrossReferences[[i]]
            
            if (is.list(ref) && !is.null(ref$database)) {
              if (ref$database == "GO") {
                go_refs <- go_refs + 1
                if (verbose && go_refs <= 3) {  # Show first 3 GO terms
                  cat("GO Reference", go_refs, ":\n")
                  cat("  ID:", ref$id, "\n")
                  if (!is.null(ref$properties)) {
                    cat("  Properties:\n")
                    for (j in seq_along(ref$properties)) {
                      prop <- ref$properties[[j]]
                      if (is.list(prop) && !is.null(prop$key) && !is.null(prop$value)) {
                        cat("    ", prop$key, ":", prop$value, "\n")
                      }
                    }
                  }
                }
              } else if (ref$database == "KEGG") {
                kegg_refs <- kegg_refs + 1
              } else {
                other_refs <- other_refs + 1
              }
            }
          }
          
          cat("Cross-reference summary:\n")
          cat("  GO terms:", go_refs, "\n")
          cat("  KEGG references:", kegg_refs, "\n")
          cat("  Other references:", other_refs, "\n")
          
        } else {
          cat("No cross-references found in API response\n")
        }
        
        # Test the funseqR extraction function
        cat("\n=== Testing funseqR extraction ===\n")
        extracted_info <- process_uniprot_json(content, debug = verbose)
        
        if (!is.null(extracted_info$go_terms)) {
          cat("funseqR extracted", nrow(extracted_info$go_terms), "GO terms\n")
          if (nrow(extracted_info$go_terms) > 0) {
            cat("First few GO terms:\n")
            print(head(extracted_info$go_terms, 3))
          }
        } else {
          cat("funseqR extracted 0 GO terms\n")
        }
        
        if (!is.null(extracted_info$kegg_refs)) {
          cat("funseqR extracted", nrow(extracted_info$kegg_refs), "KEGG references\n")
        } else {
          cat("funseqR extracted 0 KEGG references\n")
        }
        
      } else {
        cat("No results found in API response\n")
      }
      
    } else {
      cat("API call failed with status:", httr::status_code(response), "\n")
    }
    
  }, error = function(e) {
    cat("Error making API call:", e$message, "\n")
  })
}

# Test with a few problematic accessions from your results
cat("=== GO Terms Extraction Debugging ===\n")
cat("Testing specific accessions that should have GO terms...\n\n")

# Test with some known accessions from your missing GO list
test_accessions <- c("Q49GP3", "B2Z449", "E7F4Z4", "P49697")

for (acc in test_accessions) {
  debug_uniprot_response(acc, verbose = TRUE)
  cat("\n", rep("=", 80), "\n\n")
}

# Also test the evidence filter
cat("=== Testing Evidence Code Filtering ===\n")
cat("The current evidence_keep filter is:\n")
cat("c('EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC')\n")
cat("This excludes IEA (Inferred from Electronic Annotation) which is very common.\n")
cat("This might be why you're getting 0 GO terms - they might all be IEA evidence.\n\n")

# Suggest diagnostic steps
cat("=== Diagnostic Steps ===\n")
cat("1. Run this script to see the raw UniProt responses\n")
cat("2. Check if GO terms are present in the API response\n")
cat("3. Check if the funseqR extraction is working correctly\n")
cat("4. Consider relaxing the evidence code filter to include IEA\n")
cat("5. Re-run annotation with verbose=TRUE and debug_accessions\n\n")

cat("To test with relaxed evidence filtering, try:\n")
cat("annotation_results <- annotate_blast_results(con,\n")
cat("                         blast_results$blast_param_id,\n")
cat("                         evidence_keep = c('EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC', 'IEA'),\n")
cat("                         verbose = TRUE)\n")