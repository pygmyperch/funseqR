#' Extract UniProt information from API response
#'
#' @param uniprot_data The data returned from the UniProt API
#' @return A list containing extracted information or NULL if extraction fails
#' @noRd
extract_uniprot_info <- function(uniprot_data) {
  tryCatch({
    if (!is.list(uniprot_data) || is.null(uniprot_data$results) || length(uniprot_data$results) == 0) {
      warning("Unexpected or empty response from UniProt API")
      return(NULL)
    }
    
    result <- uniprot_data$results[[1]]
    
    # Extract gene names
    gene_names <- c()
    if (!is.null(result$genes)) {
      gene_names <- unlist(lapply(result$genes, function(gene) {
        c(
          if (!is.null(gene$geneName$value)) gene$geneName$value else character(0),
          if (!is.null(gene$synonyms)) sapply(gene$synonyms, function(x) x$value) else character(0),
          if (!is.null(gene$orfNames)) sapply(gene$orfNames, function(x) x$value) else character(0)
        )
      }))
      gene_names <- unique(gene_names[!is.na(gene_names)])
    }
    gene_names <- paste(gene_names, collapse = ";")
    
    # Extract GO terms
    go_terms <- data.frame(
      go_id = character(0),
      go_term = character(0),
      go_category = character(0),
      go_evidence = character(0),
      stringsAsFactors = FALSE
    )
    
    if (!is.null(result$uniProtKBCrossReferences)) {
      go_refs <- result$uniProtKBCrossReferences[sapply(result$uniProtKBCrossReferences, function(ref) ref$database == "GO")]
      if (length(go_refs) > 0) {
        go_terms <- do.call(rbind, lapply(go_refs, function(ref) {
          properties <- setNames(
            sapply(ref$properties, function(prop) prop$value),
            sapply(ref$properties, function(prop) prop$key)
          )
          data.frame(
            go_id = ref$id,
            go_term = properties["GoTerm"],
            go_category = substr(properties["GoTerm"], 1, 1),
            go_evidence = properties["GoEvidenceType"],
            stringsAsFactors = FALSE
          )
        }))
      }
    }
    
    # Extract KEGG references
    kegg_refs <- data.frame(
      kegg_id = character(0),
      stringsAsFactors = FALSE
    )
    
    if (!is.null(result$uniProtKBCrossReferences)) {
      kegg_refs <- result$uniProtKBCrossReferences[sapply(result$uniProtKBCrossReferences, function(ref) ref$database == "KEGG")]
      if (length(kegg_refs) > 0) {
        kegg_refs <- data.frame(
          kegg_id = sapply(kegg_refs, function(ref) ref$id),
          stringsAsFactors = FALSE
        )
      }
    }
    
    list(
      accession = result$primaryAccession,
      entry_name = result$uniProtkbId,
      gene_names = gene_names,
      go_terms = go_terms,
      kegg_refs = kegg_refs
    )
  }, error = function(e) {
    warning(paste("Error in extract_uniprot_info:", e$message))
    return(NULL)
  })
}
