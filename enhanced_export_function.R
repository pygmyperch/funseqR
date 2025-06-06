#!/usr/bin/env Rscript

# Enhanced export function with SNP and gene annotation details
library(funseqR)

#' Export significant GO results with detailed SNP and gene annotation information
#'
#' @param con Database connection
#' @param enrichment_results List of enrichment results from run_go_enrichment_workflow
#' @param candidate_file_id Integer. File ID of candidate dataset
#' @param background_file_id Integer. File ID of background dataset  
#' @param output_prefix Character. Prefix for output files
#' @param fdr_cutoff Numeric. FDR threshold for significance (default 0.05)
#' @param include_trending Logical. Include trending results (FDR < 0.1)
#'
export_detailed_go_results <- function(con, enrichment_results, candidate_file_id, background_file_id,
                                     output_prefix = "detailed_GO", fdr_cutoff = 0.05, 
                                     include_trending = FALSE) {
  
  # Set significance cutoff
  sig_cutoff <- if (include_trending) 0.10 else fdr_cutoff
  
  for (ontology in names(enrichment_results)) {
    if (!is.null(enrichment_results[[ontology]])) {
      
      cat("Processing", ontology, "ontology...\n")
      
      # Filter for significant terms
      sig_terms <- enrichment_results[[ontology]][
        enrichment_results[[ontology]]$p_adjusted < sig_cutoff,
      ]
      
      if (nrow(sig_terms) > 0) {
        # Sort by adjusted p-value
        sig_terms <- sig_terms[order(sig_terms$p_adjusted), ]
        
        cat("  Found", nrow(sig_terms), "significant terms\n")
        
        # For each significant GO term, get detailed SNP and gene info
        detailed_results <- list()
        
        for (i in 1:nrow(sig_terms)) {
          go_id <- sig_terms$go_id[i]
          go_term <- sig_terms$go_term[i]
          
          cat("    Processing:", go_term, "\n")
          
          # Get genes associated with this GO term in candidate dataset
          gene_query <- "
            SELECT DISTINCT
              c.chromosome,
              c.position,
              c.ref_allele,
              c.alt_allele,
              c.quality,
              c.filter,
              a.uniprot_accession,
              a.gene_names,
              a.protein_names,
              a.organism,
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
          
          gene_details <- DBI::dbGetQuery(con, gene_query, 
                                        list(candidate_file_id, background_file_id, go_id))
          
          if (nrow(gene_details) > 0) {
            # Add GO enrichment statistics to each row
            gene_details$fold_enrichment <- sig_terms$fold_enrichment[i]
            gene_details$p_value <- sig_terms$p_value[i]
            gene_details$p_adjusted <- sig_terms$p_adjusted[i]
            gene_details$significance_level <- sig_terms$significance_level[i]
            gene_details$foreground_count <- sig_terms$foreground_count[i]
            gene_details$background_count <- sig_terms$background_count[i]
            gene_details$total_foreground <- sig_terms$total_foreground[i]
            gene_details$total_background <- sig_terms$total_background[i]
            
            detailed_results[[paste0(go_id, "_", i)]] <- gene_details
          }
        }
        
        if (length(detailed_results) > 0) {
          # Combine all detailed results
          combined_details <- do.call(rbind, detailed_results)
          
          # Create summary by GO term
          go_summary <- sig_terms[, c("go_id", "go_term", "go_category", "foreground_count", 
                                    "background_count", "fold_enrichment", "p_value", 
                                    "p_adjusted", "significance_level")]
          
          # Export detailed results
          detail_filename <- paste0(output_prefix, "_", ontology, "_detailed.csv")
          write.csv(combined_details, detail_filename, row.names = FALSE)
          
          # Export summary
          summary_filename <- paste0(output_prefix, "_", ontology, "_summary.csv")
          write.csv(go_summary, summary_filename, row.names = FALSE)
          
          # Create gene-centric summary (base R approach)
          gene_keys <- paste(combined_details$uniprot_accession, 
                           combined_details$chromosome, 
                           combined_details$position, sep = "_")
          
          gene_summary_list <- list()
          for (key in unique(gene_keys)) {
            subset_data <- combined_details[gene_keys == key, ]
            
            gene_summary_list[[key]] <- data.frame(
              uniprot_accession = subset_data$uniprot_accession[1],
              gene_names = subset_data$gene_names[1],
              protein_names = subset_data$protein_names[1],
              chromosome = subset_data$chromosome[1],
              position = subset_data$position[1],
              n_significant_go_terms = nrow(subset_data),
              go_terms = paste(unique(subset_data$go_term), collapse = "; "),
              go_ids = paste(unique(subset_data$go_id), collapse = "; "),
              max_fold_enrichment = max(subset_data$fold_enrichment),
              min_p_adjusted = min(subset_data$p_adjusted),
              stringsAsFactors = FALSE
            )
          }
          
          gene_summary <- do.call(rbind, gene_summary_list)
          gene_summary <- gene_summary[order(gene_summary$min_p_adjusted), ]
          
          gene_filename <- paste0(output_prefix, "_", ontology, "_genes.csv")
          write.csv(gene_summary, gene_filename, row.names = FALSE)
          
          cat("  Exported to:\n")
          cat("    ", detail_filename, "(", nrow(combined_details), "SNP-gene-GO associations)\n")
          cat("    ", summary_filename, "(", nrow(go_summary), "GO terms)\n") 
          cat("    ", gene_filename, "(", nrow(gene_summary), "unique genes)\n")
          
        } else {
          cat("  No detailed gene information found for significant terms\n")
        }
        
      } else {
        cat("  No significant terms found (FDR <", sig_cutoff, ")\n")
      }
    }
  }
  
  cat("\nExport complete!\n")
}

# Example usage function
run_detailed_export <- function() {
  
  # Connect to database
  project_dir <- "/Users/brau0037/Library/CloudStorage/GoogleDrive-pygmyperch@gmail.com/My Drive/projects/snapperSG/south_coast/annotation/funseq_annotationMS"
  db_path <- file.path(project_dir, "funseq_project.db")
  con <- connect_funseq_db(db_path)
  
  # Run GO enrichment analysis if not already done
  cat("Running GO enrichment analysis...\n")
  results <- run_go_enrichment_workflow(
    con = con,
    project_id = 1,
    candidate_vcf_file = file.path(project_dir, "SA448_855.vcf"),
    background_file_id = 1,
    ontologies = c("BP", "MF", "CC"),
    min_genes = 3,
    store_results = TRUE,
    create_plots = FALSE,
    verbose = FALSE
  )
  
  # Export detailed results
  cat("\nExporting detailed results...\n")
  export_detailed_go_results(
    con = con,
    enrichment_results = results$enrichment_results,
    candidate_file_id = results$summary$candidate_file_id,
    background_file_id = results$summary$background_file_id,
    output_prefix = "snapper_detailed_GO",
    fdr_cutoff = 0.05,
    include_trending = TRUE  # Include trending results (FDR < 0.1)
  )
  
  close_funseq_db(con)
  cat("\nAnalysis complete! Check the generated CSV files.\n")
}

# Uncomment to run
# run_detailed_export()