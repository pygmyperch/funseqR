#!/usr/bin/env Rscript

# Export detailed GO results for snapper analysis
library(funseqR)

# Connect to database
project_dir <- "/Users/brau0037/Library/CloudStorage/GoogleDrive-pygmyperch@gmail.com/My Drive/projects/snapperSG/south_coast/annotation/funseq_annotationMS"
db_path <- file.path(project_dir, "funseq_project.db")
con <- connect_funseq_db(db_path)

cat("=== Exporting Detailed GO Results for Snapper Analysis ===\n")

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

candidate_file_id <- results$summary$candidate_file_id
background_file_id <- results$summary$background_file_id

cat("Analysis complete. Candidate genes:", results$summary$foreground_genes, "\n")
cat("Background genes:", results$summary$background_genes, "\n\n")

# Process each ontology
for (ontology in c("BP", "MF", "CC")) {
  
  if (!is.null(results$enrichment_results[[ontology]])) {
    
    cat("Processing", ontology, "ontology...\n")
    
    # Get all results (including trending)
    all_terms <- results$enrichment_results[[ontology]]
    significant_terms <- all_terms[all_terms$p_adjusted < 0.10, ]  # Include trending
    
    if (nrow(significant_terms) > 0) {
      
      cat("  Found", nrow(significant_terms), "terms with FDR < 0.10\n")
      
      # Get detailed SNP and gene information
      all_details <- list()
      
      for (i in 1:nrow(significant_terms)) {
        go_id <- significant_terms$go_id[i]
        go_term <- significant_terms$go_term[i]
        
        # Query for genes and SNPs associated with this GO term
        detail_query <- "
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
            gt.go_evidence
          FROM vcf_data c
          JOIN vcf_data r ON (c.chromosome = r.chromosome AND c.position = r.position)
          JOIN flanking_sequences fs ON r.vcf_id = fs.vcf_id
          JOIN blast_results br ON fs.flanking_id = br.flanking_id  
          JOIN annotations a ON br.blast_result_id = a.blast_result_id
          JOIN go_terms gt ON a.annotation_id = gt.annotation_id
          WHERE c.file_id = ? AND r.file_id = ? AND gt.go_id = ?
          ORDER BY c.chromosome, c.position
        "
        
        gene_details <- DBI::dbGetQuery(con, detail_query, 
                                      list(candidate_file_id, background_file_id, go_id))
        
        if (nrow(gene_details) > 0) {
          # Add enrichment statistics
          gene_details$fold_enrichment <- significant_terms$fold_enrichment[i]
          gene_details$p_value <- significant_terms$p_value[i]
          gene_details$p_adjusted <- significant_terms$p_adjusted[i]
          gene_details$significance_level <- significant_terms$significance_level[i]
          gene_details$foreground_count <- significant_terms$foreground_count[i]
          gene_details$background_count <- significant_terms$background_count[i]
          
          all_details[[go_id]] <- gene_details
        }
      }
      
      if (length(all_details) > 0) {
        # Combine all results
        combined_data <- do.call(rbind, all_details)
        
        # Export main detailed file
        detail_file <- paste0("snapper_", ontology, "_detailed_results.csv")
        write.csv(combined_data, detail_file, row.names = FALSE)
        
        # Create GO term summary
        go_summary <- significant_terms[order(significant_terms$p_adjusted), ]
        summary_file <- paste0("snapper_", ontology, "_GO_summary.csv")
        write.csv(go_summary, summary_file, row.names = FALSE)
        
        # Create gene summary
        gene_keys <- paste(combined_data$chromosome, combined_data$position, 
                          combined_data$uniprot_accession, sep = "_")
        
        gene_list <- list()
        for (key in unique(gene_keys)) {
          subset_data <- combined_data[gene_keys == key, ]
          
          gene_list[[key]] <- data.frame(
            chromosome = subset_data$chromosome[1],
            position = subset_data$position[1],
            ref_allele = subset_data$ref[1],
            alt_allele = subset_data$alt[1],
            uniprot_accession = subset_data$uniprot_accession[1],
            gene_names = subset_data$gene_names[1],
            entry_name = subset_data$entry_name[1],
            n_enriched_go_terms = length(unique(subset_data$go_id)),
            enriched_go_terms = paste(unique(subset_data$go_term), collapse = " | "),
            max_fold_enrichment = max(subset_data$fold_enrichment),
            best_p_adjusted = min(subset_data$p_adjusted),
            best_e_value = min(subset_data$e_value),
            max_bit_score = max(subset_data$bit_score),
            stringsAsFactors = FALSE
          )
        }
        
        gene_summary <- do.call(rbind, gene_list)
        gene_summary <- gene_summary[order(gene_summary$best_p_adjusted), ]
        gene_file <- paste0("snapper_", ontology, "_gene_summary.csv")
        write.csv(gene_summary, gene_file, row.names = FALSE)
        
        # Print summary
        highly_sig <- sum(significant_terms$p_adjusted < 0.01)
        sig <- sum(significant_terms$p_adjusted < 0.05)
        trending <- sum(significant_terms$p_adjusted < 0.10) - sig
        
        cat("  Results exported:\n")
        cat("    ", detail_file, "-", nrow(combined_data), "SNP-gene-GO associations\n")
        cat("    ", summary_file, "-", nrow(go_summary), "GO terms\n")
        cat("    ", gene_file, "-", nrow(gene_summary), "unique genes\n")
        cat("  Significance breakdown:\n")
        cat("    Highly significant (FDR < 0.01):", highly_sig, "terms\n")
        cat("    Significant (FDR < 0.05):", sig, "terms\n")
        cat("    Trending (FDR < 0.10):", trending, "terms\n\n")
        
      } else {
        cat("  No gene details found for significant terms\n\n")
      }
      
    } else {
      cat("  No significant terms found (FDR < 0.10)\n\n")
    }
    
  } else {
    cat("No results for", ontology, "ontology\n\n")
  }
}

close_funseq_db(con)

cat("=== Export Complete ===\n")
cat("Files generated:\n")
files <- list.files(pattern = "snapper_.+\\.csv$")
for (file in files) {
  cat("  ", file, "\n")
}
cat("\nThese files contain:\n")
cat("  *_detailed_results.csv: Complete SNP, gene, and GO term associations\n")
cat("  *_GO_summary.csv: GO term enrichment statistics\n")
cat("  *_gene_summary.csv: Gene-centric view with enrichment summaries\n")