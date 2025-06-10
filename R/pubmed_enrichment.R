#' PubMed Literature Search for Enriched Genes
#'
#' Functions to search PubMed for literature related to genes associated with
#' enriched GO terms and generate reference lists.

#' Enhanced PubMed search for genes with multiple strategies
#'
#' @param uniprot_accession Character. UniProt accession number
#' @param gene_name Character. Gene name or symbol (can be empty)
#' @param entry_name Character. UniProt entry name
#' @param species Character. Species name (default "danio rerio")
#' @param max_results Integer. Maximum number of papers to retrieve (default 5)
#' @param custom_search_term Character. Optional custom search term to include
#' @param include_functional_mapping Logical. Use known functional mappings (default TRUE)
#' @param verbose Logical. Print progress information (default TRUE)
#'
#' @return Data frame with PubMed search results
#'
#' @details
#' This function uses multiple search strategies to find relevant literature:
#' 1. Gene names (if available and suitable length)
#' 2. UniProt accessions
#' 3. Functional names derived from entry names
#' 4. Custom search terms provided by user
#' 5. Known functional mappings (lectin, chymotrypsinogen, etc.)
#'
#' @export
search_pubmed_enhanced <- function(uniprot_accession, gene_name = "", entry_name = "",
                                 species = "danio rerio", max_results = 5,
                                 custom_search_term = NULL, include_functional_mapping = TRUE,
                                 verbose = TRUE) {
  
  if (verbose) message("Searching PubMed for: ", uniprot_accession, 
                      if(gene_name != "") paste0(" (", gene_name, ")") else "")
  
  # Create multiple search strategies
  search_terms <- c()
  
  # Strategy 1: Use gene name if available and reasonable length
  if (!is.na(gene_name) && gene_name != "" && nchar(gene_name) >= 3 && nchar(gene_name) <= 10) {
    search_terms <- c(search_terms, gene_name)
  }
  
  # Strategy 2: Use UniProt accession
  if (!is.na(uniprot_accession) && uniprot_accession != "") {
    search_terms <- c(search_terms, uniprot_accession)
  }
  
  # Strategy 3: Extract meaningful parts from entry name
  if (!is.na(entry_name) && entry_name != "") {
    # Extract the prefix before underscore (often the functional name)
    entry_parts <- strsplit(entry_name, "_")[[1]]
    if (length(entry_parts) > 0 && nchar(entry_parts[1]) >= 3) {
      search_terms <- c(search_terms, entry_parts[1])
    }
  }
  
  # Strategy 4: Custom search term provided by user
  if (!is.null(custom_search_term) && custom_search_term != "") {
    search_terms <- c(search_terms, custom_search_term)
  }
  
  # Strategy 5: Known functional mappings
  if (include_functional_mapping) {
    functional_terms <- list(
      "LEC" = "lectin",
      "CTRB" = "chymotrypsinogen B",
      "CTRA" = "chymotrypsinogen A", 
      "CTXB" = "cytotoxin",
      "STXA" = "stonustoxin",
      "STXB" = "stonustoxin",
      "FUCL" = "fucolectin",
      "GLRA" = "glycine receptor",
      "ASAH" = "acid ceramidase",
      "ELAVL" = "ELAV like",
      "KLHL" = "kelch like",
      "TENM" = "teneurin",
      "MYO1" = "myosin",
      "POC1" = "POC1 centriolar protein"
    )
    
    for (term in names(functional_terms)) {
      if (grepl(term, entry_name, ignore.case = TRUE)) {
        search_terms <- c(search_terms, functional_terms[[term]])
      }
    }
  }
  
  # Remove duplicates and empty terms
  search_terms <- unique(search_terms[search_terms != "" & !is.na(search_terms)])
  
  if (length(search_terms) == 0) {
    if (verbose) message("  No valid search terms found")
    return(data.frame())
  }
  
  if (verbose) message("  Search terms: ", paste(search_terms, collapse = ", "))
  
  all_papers <- list()
  
  for (term in search_terms) {
    
    # Create species-specific search terms
    species_terms <- c()
    if (species != "") {
      species_parts <- strsplit(species, " ")[[1]]
      if (length(species_parts) >= 2) {
        genus <- species_parts[1]
        species_name <- species_parts[2]
        species_terms <- c(
          paste0('"', species, '"'),
          genus,
          if (genus == "danio") c("zebrafish", "fish") else c("fish"),
          if (!is.null(custom_search_term)) custom_search_term else NULL
        )
      }
    }
    
    # Try multiple search strategies for each term
    search_queries <- c()
    
    # Query A: Specific with full species name
    if (length(species_terms) > 0) {
      search_queries <- c(search_queries,
        paste0('"', term, '"[Title/Abstract] AND ', species_terms[1], '[MeSH Terms OR Title/Abstract]')
      )
    }
    
    # Query B: Broader with fish terms
    search_queries <- c(search_queries,
      paste0('"', term, '"[Title/Abstract] AND ("fish"[Title/Abstract] OR "teleost"[Title/Abstract])')
    )
    
    # Query C: UniProt accession specific search (if term is the accession)
    if (term == uniprot_accession) {
      search_queries <- c(search_queries,
        paste0('"', term, '"[All Fields]')
      )
    }
    
    # Query D: Very broad functional search
    search_queries <- c(search_queries,
      paste0('"', term, '"[Title/Abstract]')
    )
    
    # Query E: Custom search term combination
    if (!is.null(custom_search_term) && custom_search_term != "" && term != custom_search_term) {
      search_queries <- c(search_queries,
        paste0('"', term, '"[Title/Abstract] AND "', custom_search_term, '"[Title/Abstract]')
      )
    }
    
    for (query in search_queries) {
      
      tryCatch({
        
        # URL encode the search term
        encoded_search <- URLencode(query, reserved = TRUE)
        
        # NCBI E-utilities URLs
        esearch_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?",
                             "db=pubmed&term=", encoded_search, 
                             "&retmax=", max_results, "&retmode=xml")
        
        # Search for PMIDs
        search_response <- xml2::read_xml(esearch_url)
        pmids <- xml2::xml_text(xml2::xml_find_all(search_response, "//Id"))
        
        if (length(pmids) > 0) {
          if (verbose) message("    Found ", length(pmids), " papers with: ", substr(query, 1, 60), "...")
          
          # Get paper details
          pmid_list <- paste(pmids, collapse = ",")
          efetch_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?",
                              "db=pubmed&id=", pmid_list, "&retmode=xml")
          
          Sys.sleep(0.5)  # Be respectful to NCBI servers
          
          fetch_response <- xml2::read_xml(efetch_url)
          articles <- xml2::xml_find_all(fetch_response, "//PubmedArticle")
          
          for (i in seq_along(articles)) {
            article <- articles[[i]]
            
            pmid <- xml2::xml_text(xml2::xml_find_first(article, ".//PMID"))
            title <- xml2::xml_text(xml2::xml_find_first(article, ".//ArticleTitle"))
            
            # Authors
            author_nodes <- xml2::xml_find_all(article, ".//Author")
            authors <- sapply(author_nodes, function(auth) {
              last <- xml2::xml_text(xml2::xml_find_first(auth, ".//LastName"))
              first <- xml2::xml_text(xml2::xml_find_first(auth, ".//ForeName"))
              if (!is.na(last) && !is.na(first)) {
                paste(last, substr(first, 1, 1), sep = " ")
              } else if (!is.na(last)) {
                last
              } else {
                NA
              }
            })
            authors <- authors[!is.na(authors)]
            author_string <- if (length(authors) > 0) {
              if (length(authors) > 3) {
                paste(paste(authors[1:3], collapse = ", "), "et al.")
              } else {
                paste(authors, collapse = ", ")
              }
            } else {
              "No authors listed"
            }
            
            # Journal and year
            journal <- xml2::xml_text(xml2::xml_find_first(article, ".//Journal/Title"))
            year <- xml2::xml_text(xml2::xml_find_first(article, ".//PubDate/Year"))
            
            # Abstract (first 300 characters)
            abstract <- xml2::xml_text(xml2::xml_find_first(article, ".//Abstract/AbstractText"))
            if (!is.na(abstract) && nchar(abstract) > 300) {
              abstract <- paste0(substr(abstract, 1, 300), "...")
            }
            
            paper_info <- data.frame(
              uniprot_accession = uniprot_accession,
              gene_name = gene_name,
              entry_name = entry_name,
              search_term = term,
              search_query = substr(query, 1, 100),
              pmid = pmid %||% NA,
              title = title %||% NA,
              authors = author_string,
              journal = journal %||% NA,
              year = year %||% NA,
              abstract = abstract %||% NA,
              pubmed_url = paste0("https://pubmed.ncbi.nlm.nih.gov/", pmid %||% ""),
              stringsAsFactors = FALSE
            )
            
            all_papers[[paste0(pmid, "_", term)]] <- paper_info
          }
          
          # If we found papers with this query, don't try more aggressive searches for this term
          if (length(pmids) >= 2) break
        }
        
        Sys.sleep(0.3)  # Be respectful
        
      }, error = function(e) {
        if (verbose) message("    Error with query: ", e$message)
      })
    }
    
    # If we found papers with this term, don't try other terms unless we have very few
    if (length(all_papers) >= 3) break
  }
  
  if (length(all_papers) > 0) {
    result <- do.call(rbind, all_papers)
    # Remove duplicates by PMID
    result <- result[!duplicated(result$pmid), ]
    
    # Sort by year (newest first)
    if (nrow(result) > 0 && !all(is.na(result$year))) {
      result$year_numeric <- as.numeric(result$year)
      result <- result[order(-result$year_numeric, na.last = TRUE), ]
      result$year_numeric <- NULL
    }
    
    return(result)
  } else {
    if (verbose) message("  No papers found")
    return(data.frame())
  }
}

# Helper function for null coalescing
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x

#' Get literature for enriched genes with enhanced search
#'
#' @param con Database connection
#' @param candidate_file_id Integer. File ID of candidate dataset
#' @param background_file_id Integer. File ID of background dataset
#' @param enriched_go_ids Vector. GO IDs of enriched terms
#' @param species Character. Species name for PubMed search (default "danio rerio")
#' @param max_papers_per_gene Integer. Maximum papers per gene (default 3)
#' @param custom_search_term Character. Optional custom search term to include with all searches
#' @param include_functional_mapping Logical. Use known functional mappings (default TRUE)
#' @param verbose Logical. Print progress information (default TRUE)
#'
#' @return List containing gene information and literature references
#'
#' @export
get_literature_for_enriched_genes <- function(con, candidate_file_id, background_file_id,
                                             enriched_go_ids, species = "danio rerio",
                                             max_papers_per_gene = 3, custom_search_term = NULL,
                                             include_functional_mapping = TRUE, verbose = TRUE) {
  
  if (verbose) message("Getting literature for enriched genes...")
  
  # Get unique genes associated with enriched GO terms
  gene_query <- paste(
    "SELECT",
    "  a.uniprot_accession,",
    "  a.gene_names,",
    "  a.entry_name,",
    "  COUNT(DISTINCT gt.go_id) as n_enriched_go_terms,",
    "  GROUP_CONCAT(gt.go_term, ' | ') as go_terms",
    "FROM vcf_data c",
    "JOIN vcf_data r ON (c.chromosome = r.chromosome AND c.position = r.position)",
    "JOIN flanking_sequences fs ON r.vcf_id = fs.vcf_id",
    "JOIN blast_results br ON fs.flanking_id = br.flanking_id",
    "JOIN annotations a ON br.blast_result_id = a.blast_result_id",
    "JOIN go_terms gt ON a.annotation_id = gt.annotation_id",
    "WHERE c.file_id = ? AND r.file_id = ?",
    "  AND gt.go_id IN (", paste0("'", enriched_go_ids, "'", collapse = ", "), ")",
    "GROUP BY a.uniprot_accession, a.gene_names, a.entry_name",
    "ORDER BY n_enriched_go_terms DESC, a.gene_names"
  )
  
  genes <- DBI::dbGetQuery(con, gene_query, list(candidate_file_id, background_file_id))
  
  if (nrow(genes) == 0) {
    if (verbose) message("No genes found for enriched GO terms")
    return(list(genes = data.frame(), literature = data.frame()))
  }
  
  if (verbose) message("Found ", nrow(genes), " unique genes")
  
  if (verbose) message("Searching literature for ", nrow(genes), " genes")
  
  # Search PubMed for each gene using enhanced search
  all_literature <- list()
  
  for (i in 1:nrow(genes)) {
    
    if (verbose) message("  [", i, "/", nrow(genes), "] ", genes$uniprot_accession[i])
    
    papers <- search_pubmed_enhanced(
      uniprot_accession = genes$uniprot_accession[i],
      gene_name = genes$gene_names[i],
      entry_name = genes$entry_name[i],
      species = species,
      max_results = max_papers_per_gene,
      custom_search_term = custom_search_term,
      include_functional_mapping = include_functional_mapping,
      verbose = FALSE
    )
    
    if (nrow(papers) > 0) {
      # Add gene metadata from the original query
      if ("n_enriched_go_terms" %in% names(genes)) {
        papers$n_enriched_go_terms <- genes$n_enriched_go_terms[i]
      }
      if ("go_terms" %in% names(genes)) {
        papers$associated_go_terms <- genes$go_terms[i]
      }
      
      all_literature[[genes$uniprot_accession[i]]] <- papers
    }
    
    # Be respectful to NCBI servers
    Sys.sleep(1)
  }
  
  # Combine literature results
  if (length(all_literature) > 0) {
    literature_df <- do.call(rbind, all_literature)
    rownames(literature_df) <- NULL
  } else {
    literature_df <- data.frame()
  }
  
  if (verbose) {
    message("Literature search complete:")
    message("  Genes searched: ", nrow(genes))
    message("  Papers found: ", nrow(literature_df))
    message("  Genes with literature: ", length(unique(literature_df$uniprot_accession)))
  }
  
  return(list(
    genes = genes,
    literature = literature_df,
    summary = list(
      total_genes = nrow(genes),
      papers_found = nrow(literature_df),
      genes_with_literature = length(unique(literature_df$uniprot_accession))
    )
  ))
}

#' Generate formatted reference list
#'
#' @param literature_data Data frame. Output from get_literature_for_enriched_genes()$literature
#' @param output_file Character. Optional file to write references
#' @param format Character. Reference format: "apa", "basic", or "detailed"
#'
#' @return Character vector of formatted references
#'
#' @export
generate_reference_list <- function(literature_data, output_file = NULL, format = "apa") {
  
  if (nrow(literature_data) == 0) {
    message("No literature data to format")
    return(character(0))
  }
  
  references <- c()
  
  # Group by gene
  for (gene in unique(literature_data$gene)) {
    gene_papers <- literature_data[literature_data$gene == gene, ]
    
    references <- c(references, paste("\n## References for", gene, "\n"))
    
    for (i in 1:nrow(gene_papers)) {
      paper <- gene_papers[i, ]
      
      if (format == "apa") {
        # APA format
        ref <- paste0(
          paper$authors, " (", paper$year, "). ",
          paper$title, ". ",
          paper$journal, ". ",
          "PMID: ", paper$pmid, ". ",
          "https://pubmed.ncbi.nlm.nih.gov/", paper$pmid
        )
      } else if (format == "basic") {
        # Basic format
        ref <- paste0(
          paper$title, " (", paper$year, "). ",
          paper$authors, ". ", paper$journal, ". ",
          "PMID: ", paper$pmid
        )
      } else {
        # Detailed format
        ref <- paste0(
          "**", paper$title, "**\n",
          "Authors: ", paper$authors, "\n",
          "Journal: ", paper$journal, " (", paper$year, ")\n",
          "PMID: ", paper$pmid, "\n",
          "URL: https://pubmed.ncbi.nlm.nih.gov/", paper$pmid, "\n",
          if (!is.na(paper$abstract)) paste0("Abstract: ", paper$abstract, "\n"),
          "Associated GO terms: ", paper$associated_go_terms, "\n"
        )
      }
      
      references <- c(references, paste(i, ".", ref, "\n"))
    }
  }
  
  # Write to file if requested
  if (!is.null(output_file)) {
    writeLines(references, output_file)
    message("References written to: ", output_file)
  }
  
  return(references)
}