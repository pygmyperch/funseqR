#!/usr/bin/env Rscript

# Export complete functional annotations (GO + KEGG) for snapper analysis
library(funseqR)

# Connect to database
project_dir <- "/Users/brau0037/Library/CloudStorage/GoogleDrive-pygmyperch@gmail.com/My Drive/projects/snapperSG/south_coast/annotation/funseq_annotationMS"
db_path <- file.path(project_dir, "funseq_project.db")
con <- connect_funseq_db(db_path)

cat("=== Exporting Complete Functional Annotations (GO + KEGG) ===\n")

# First, let's get the complete functional annotation picture
cat("Getting complete functional annotation data...\n")

# Query for all candidate SNPs with their complete functional annotations
complete_query <- "
  SELECT DISTINCT
    c.chromosome,
    c.position,
    c.ref,
    c.alt,
    c.qual,
    c.filter,
    a.uniprot_accession,
    a.gene_names,
    a.entry_name,
    br.e_value,
    br.bit_score,
    br.percent_identity,
    gt.go_id,
    gt.go_term,
    gt.go_category,
    gt.go_evidence,
    kr.kegg_id,
    kr.pathway_name
  FROM vcf_data c
  JOIN vcf_data r ON (c.chromosome = r.chromosome AND c.position = r.position)
  JOIN flanking_sequences fs ON r.vcf_id = fs.vcf_id
  JOIN blast_results br ON fs.flanking_id = br.flanking_id  
  JOIN annotations a ON br.blast_result_id = a.blast_result_id
  LEFT JOIN go_terms gt ON a.annotation_id = gt.annotation_id
  LEFT JOIN kegg_references kr ON a.annotation_id = kr.annotation_id
  WHERE c.file_id = ? AND r.file_id = ?
  ORDER BY c.chromosome, c.position, a.uniprot_accession
"

candidate_file_id <- 3  # SA448_855.vcf
background_file_id <- 1  # SA448_14699.vcf

cat("Extracting complete annotation data...\n")
complete_annotations <- DBI::dbGetQuery(con, complete_query, 
                                      list(candidate_file_id, background_file_id))

cat("Found", nrow(complete_annotations), "annotation records\n")

# Run GO enrichment analysis
cat("Running GO enrichment analysis...\n")
results <- run_go_enrichment_workflow(
  con = con,
  project_id = 1,
  candidate_vcf_file = file.path(project_dir, "SA448_855.vcf"),
  background_file_id = 1,
  ontologies = c("BP", "MF", "CC"),
  min_genes = 3,
  store_results = FALSE,
  create_plots = FALSE,
  verbose = FALSE
)

cat("GO analysis complete.\n")
cat("Candidate genes:", results$summary$foreground_genes, "\n")
cat("Background genes:", results$summary$background_genes, "\n")

# Export 1: Complete functional annotation table
cat("\n1. Exporting complete functional annotation table...\n")
complete_file <- "snapper_complete_functional_annotations.csv"
write.csv(complete_annotations, complete_file, row.names = FALSE)

# Create gene-centric summary with all functional info
cat("2. Creating gene-centric functional summary...\n")

# Group by gene and position to summarize all functional annotations
gene_keys <- paste(complete_annotations$chromosome, 
                  complete_annotations$position, 
                  complete_annotations$uniprot_accession, sep = "_")

gene_functional_summary <- list()
for (key in unique(gene_keys)) {
  subset_data <- complete_annotations[gene_keys == key, ]
  
  # Get unique GO terms and KEGG pathways
  go_terms <- subset_data[!is.na(subset_data$go_term) & subset_data$go_term != "", ]
  kegg_terms <- subset_data[!is.na(subset_data$kegg_id) & subset_data$kegg_id != "", ]
  
  # Separate GO by category
  bp_terms <- go_terms[go_terms$go_category == "P", ]
  mf_terms <- go_terms[go_terms$go_category == "F", ]
  cc_terms <- go_terms[go_terms$go_category == "C", ]
  
  gene_functional_summary[[key]] <- data.frame(
    chromosome = subset_data$chromosome[1],
    position = subset_data$position[1],
    ref_allele = subset_data$ref[1],
    alt_allele = subset_data$alt[1],
    quality = subset_data$qual[1],
    uniprot_accession = subset_data$uniprot_accession[1],
    gene_names = subset_data$gene_names[1],
    entry_name = subset_data$entry_name[1],
    best_e_value = min(subset_data$e_value, na.rm = TRUE),
    max_bit_score = max(subset_data$bit_score, na.rm = TRUE),
    max_percent_identity = max(subset_data$percent_identity, na.rm = TRUE),
    
    # GO term counts
    n_go_terms_total = nrow(go_terms),
    n_bp_terms = nrow(bp_terms),
    n_mf_terms = nrow(mf_terms),
    n_cc_terms = nrow(cc_terms),
    
    # GO terms (separated by category)
    biological_processes = if(nrow(bp_terms) > 0) paste(unique(bp_terms$go_term), collapse = " | ") else "",
    molecular_functions = if(nrow(mf_terms) > 0) paste(unique(mf_terms$go_term), collapse = " | ") else "",
    cellular_components = if(nrow(cc_terms) > 0) paste(unique(cc_terms$go_term), collapse = " | ") else "",
    
    # KEGG information
    n_kegg_pathways = length(unique(kegg_terms$kegg_id[kegg_terms$kegg_id != ""])),
    kegg_ids = if(nrow(kegg_terms) > 0) paste(unique(kegg_terms$kegg_id), collapse = " | ") else "",
    kegg_pathways = if(nrow(kegg_terms) > 0) paste(unique(kegg_terms$pathway_name[kegg_terms$pathway_name != "-"]), collapse = " | ") else "",
    
    stringsAsFactors = FALSE
  )
}

gene_summary <- do.call(rbind, gene_functional_summary)
gene_summary_file <- "snapper_gene_functional_summary.csv"
write.csv(gene_summary, gene_summary_file, row.names = FALSE)

# Export 3: GO enrichment results with associated genes and KEGG info
cat("3. Exporting GO enrichment results with KEGG cross-references...\n")

for (ontology in c("BP", "MF", "CC")) {
  
  if (!is.null(results$enrichment_results[[ontology]])) {
    
    cat("  Processing", ontology, "ontology...\n")
    
    # Get significant/trending GO terms
    go_results <- results$enrichment_results[[ontology]]
    sig_go <- go_results[go_results$p_adjusted < 0.10, ]  # Include trending
    
    if (nrow(sig_go) > 0) {
      
      # For each significant GO term, get associated genes and their KEGG info
      enriched_with_kegg <- list()
      
      for (i in 1:nrow(sig_go)) {
        go_id <- sig_go$go_id[i]
        
        # Get genes associated with this GO term
        go_gene_query <- "
          SELECT DISTINCT
            c.chromosome,
            c.position,
            c.ref,
            c.alt,
            a.uniprot_accession,
            a.gene_names,
            a.entry_name,
            br.e_value,
            br.bit_score,
            br.percent_identity,
            gt.go_id,
            gt.go_term,
            gt.go_category,
            kr.kegg_id,
            kr.pathway_name
          FROM vcf_data c
          JOIN vcf_data r ON (c.chromosome = r.chromosome AND c.position = r.position)
          JOIN flanking_sequences fs ON r.vcf_id = fs.vcf_id
          JOIN blast_results br ON fs.flanking_id = br.flanking_id  
          JOIN annotations a ON br.blast_result_id = a.blast_result_id
          JOIN go_terms gt ON a.annotation_id = gt.annotation_id
          LEFT JOIN kegg_references kr ON a.annotation_id = kr.annotation_id
          WHERE c.file_id = ? AND r.file_id = ? AND gt.go_id = ?
          ORDER BY c.chromosome, c.position
        "
        
        gene_data <- DBI::dbGetQuery(con, go_gene_query, 
                                   list(candidate_file_id, background_file_id, go_id))
        
        if (nrow(gene_data) > 0) {
          # Add GO enrichment statistics
          gene_data$fold_enrichment <- sig_go$fold_enrichment[i]
          gene_data$p_value <- sig_go$p_value[i]
          gene_data$p_adjusted <- sig_go$p_adjusted[i]
          gene_data$significance_level <- sig_go$significance_level[i]
          gene_data$foreground_count <- sig_go$foreground_count[i]
          gene_data$background_count <- sig_go$background_count[i]
          
          enriched_with_kegg[[go_id]] <- gene_data
        }
      }
      
      if (length(enriched_with_kegg) > 0) {
        # Combine all enriched results
        combined_enriched <- do.call(rbind, enriched_with_kegg)
        
        # Export the detailed enrichment file
        enriched_file <- paste0("snapper_", ontology, "_enriched_with_kegg.csv")
        write.csv(combined_enriched, enriched_file, row.names = FALSE)
        
        cat("    Exported", enriched_file, "with", nrow(combined_enriched), "records\n")
      }
    } else {
      cat("    No significant", ontology, "terms found\n")
    }
  }
}

# Export 4: KEGG pathway summary
cat("4. Creating KEGG pathway summary...\n")

# Get KEGG pathway frequencies in candidate genes
kegg_summary_query <- "
  SELECT 
    kr.kegg_id,
    kr.pathway_name,
    COUNT(DISTINCT c.chromosome || '_' || c.position) as n_snps,
    COUNT(DISTINCT a.uniprot_accession) as n_genes,
    GROUP_CONCAT(DISTINCT a.gene_names) as gene_list
  FROM vcf_data c
  JOIN vcf_data r ON (c.chromosome = r.chromosome AND c.position = r.position)
  JOIN flanking_sequences fs ON r.vcf_id = fs.vcf_id
  JOIN blast_results br ON fs.flanking_id = br.flanking_id  
  JOIN annotations a ON br.blast_result_id = a.blast_result_id
  JOIN kegg_references kr ON a.annotation_id = kr.annotation_id
  WHERE c.file_id = ? AND r.file_id = ? 
    AND kr.kegg_id != '' AND kr.kegg_id IS NOT NULL
  GROUP BY kr.kegg_id, kr.pathway_name
  ORDER BY n_genes DESC, n_snps DESC
"

kegg_summary <- DBI::dbGetQuery(con, kegg_summary_query, 
                               list(candidate_file_id, background_file_id))

if (nrow(kegg_summary) > 0) {
  kegg_file <- "snapper_kegg_pathway_summary.csv"
  write.csv(kegg_summary, kegg_file, row.names = FALSE)
  cat("  Exported", kegg_file, "with", nrow(kegg_summary), "KEGG pathways\n")
} else {
  cat("  No KEGG pathway data found\n")
}

close_funseq_db(con)

cat("\n=== Export Complete ===\n")
cat("Files generated:\n")
files <- list.files(pattern = "snapper_.+\\.csv$")
for (file in files) {
  file_info <- file.info(file)
  cat("  ", file, "(", round(file_info$size/1024, 1), "KB )\n")
}

cat("\nFile descriptions:\n")
cat("  snapper_complete_functional_annotations.csv: All SNP-gene-GO-KEGG associations\n")
cat("  snapper_gene_functional_summary.csv: Gene-centric view with all functional annotations\n")
cat("  snapper_*_enriched_with_kegg.csv: GO enriched terms with associated KEGG pathways\n")
cat("  snapper_kegg_pathway_summary.csv: KEGG pathway frequencies in candidate genes\n")