summarize_annotations <- function(annotations) {
  library(dplyr)
  library(tidyr)
  library(stringr)
  
  # Start with the annotated BLAST results
  summary <- annotations$annotated_blast %>%
    dplyr::select(qseqid, sseqid, uniprot_accession, gene_names, pident, length, evalue, bitscore)
  
  # Add KEGG information if available
  if (!is.null(annotations$kegg_refs)) {
    kegg_info <- annotations$kegg_refs %>%
      dplyr::select(qseqid, kegg_id) %>%
      distinct() %>%
      group_by(qseqid) %>%
      summarise(kegg_id = paste(unique(kegg_id), collapse = "; "))
    
    summary <- summary %>%
      left_join(kegg_info, by = "qseqid")
  } else {
    summary$kegg_id <- NA_character_
  }
  
  # Add GO terms if available
  if (!is.null(annotations$go_terms)) {
    go_summary <- annotations$go_terms %>%
      group_by(qseqid) %>%
      summarise(
        go_ids = paste(unique(go_id), collapse = "; "),
        go_terms = paste(unique(go_term), collapse = "; "),
        go_categories = paste(sort(unique(go_category)), collapse = "; ")
      )
    
    summary <- summary %>%
      left_join(go_summary, by = "qseqid")
  } else {
    summary$go_ids <- NA_character_
    summary$go_terms <- NA_character_
    summary$go_categories <- NA_character_
  }
  
  # Extract protein names from sseqid
  summary <- summary %>%
    mutate(
      protein_name = str_extract(sseqid, "(?<=\\|)[^_]+(?=_)")
    )
  
  # Rename columns for clarity
  summary <- summary %>%
    rename(
      snp_id = qseqid,
      uniprot_entry = sseqid
    )
  
  # Reorder columns
  summary <- summary %>%
    dplyr::select(snp_id, uniprot_entry, uniprot_accession, kegg_id, gene_names, protein_name, pident, length, evalue, bitscore, go_ids, go_terms, go_categories, everything())
  
  return(summary)
}