#!/usr/bin/env Rscript

# Comprehensive ORA Gene Set Composition Analysis
# Comparing cp_enrich_funseq_project.db vs simple_funseq_project.db

library(funseqR)
library(DBI)

cat("=== COMPREHENSIVE ORA ANALYSIS: Gene Set Composition Investigation ===\n")
cat("Comparing cp_enrich_funseq_project.db vs simple_funseq_project.db\n\n")

# Connect to both databases
cat("1. CONNECTING TO DATABASES\n")
cat("----------------------------\n")

con_cp <- tryCatch({
  connect_funseq_db('cp_enrich_funseq_project.db')
}, error = function(e) {
  cat("ERROR connecting to cp_enrich database:", e$message, "\n")
  return(NULL)
})

con_simple <- tryCatch({
  connect_funseq_db('simple_funseq_project.db')
}, error = function(e) {
  cat("ERROR connecting to simple database:", e$message, "\n")
  return(NULL)
})

if (is.null(con_cp) || is.null(con_simple)) {
  cat("Cannot proceed without both database connections\n")
  quit(status = 1)
}

cat("✓ Both databases connected successfully\n\n")

# Function to get database stats
get_db_stats <- function(con, db_name) {
  cat("Database:", db_name, "\n")
  
  # Basic table counts
  tables <- DBI::dbListTables(con)
  cat("  Tables available:", paste(tables, collapse = ", "), "\n")
  
  # Total annotations
  total_annotations <- DBI::dbGetQuery(con, "SELECT COUNT(*) as count FROM annotations")$count
  cat("  Total annotations:", total_annotations, "\n")
  
  # Annotations with GO terms
  go_annotations <- DBI::dbGetQuery(con, "
    SELECT COUNT(DISTINCT a.annotation_id) as count 
    FROM annotations a 
    JOIN go_terms gt ON a.annotation_id = gt.annotation_id
  ")$count
  cat("  Annotations with GO terms:", go_annotations, "\n")
  
  # Annotations with KEGG terms
  kegg_annotations <- DBI::dbGetQuery(con, "
    SELECT COUNT(DISTINCT a.annotation_id) as count 
    FROM annotations a 
    JOIN kegg_references kr ON a.annotation_id = kr.annotation_id
  ")$count
  cat("  Annotations with KEGG terms:", kegg_annotations, "\n")
  
  # Candidate loci count
  if ("candidate_loci" %in% tables) {
    candidate_count <- DBI::dbGetQuery(con, "SELECT COUNT(*) as count FROM candidate_loci")$count
    cat("  Candidate loci:", candidate_count, "\n")
  } else {
    cat("  Candidate loci: NO TABLE\n")
  }
  
  # VCF data count
  if ("vcf_data" %in% tables) {
    vcf_count <- DBI::dbGetQuery(con, "SELECT COUNT(*) as count FROM vcf_data")$count
    cat("  VCF data entries:", vcf_count, "\n")
  }
  
  cat("\n")
  
  return(list(
    total_annotations = total_annotations,
    go_annotations = go_annotations,
    kegg_annotations = kegg_annotations,
    candidate_count = if ("candidate_loci" %in% tables) candidate_count else 0
  ))
}

# Get basic statistics
cat("2. DATABASE BASIC STATISTICS\n")
cat("==============================\n")

stats_cp <- get_db_stats(con_cp, "cp_enrich_funseq_project.db")
stats_simple <- get_db_stats(con_simple, "simple_funseq_project.db")

# Compare candidate loci
cat("3. CANDIDATE LOCI COMPARISON\n")
cat("==============================\n")

# Get candidate loci from both databases
candidates_cp <- tryCatch({
  DBI::dbGetQuery(con_cp, "
    SELECT chromosome, position, 
           chromosome || '_' || position as locus_key
    FROM candidate_loci 
    ORDER BY chromosome, position
  ")
}, error = function(e) {
  cat("Error getting cp_enrich candidates:", e$message, "\n")
  data.frame()
})

candidates_simple <- tryCatch({
  DBI::dbGetQuery(con_simple, "
    SELECT chromosome, position,
           chromosome || '_' || position as locus_key  
    FROM candidate_loci 
    ORDER BY chromosome, position
  ")
}, error = function(e) {
  cat("Error getting simple candidates:", e$message, "\n")
  data.frame()
})

cat("cp_enrich candidates:", nrow(candidates_cp), "\n")
cat("simple candidates:", nrow(candidates_simple), "\n")

if (nrow(candidates_cp) > 0 && nrow(candidates_simple) > 0) {
  # Check if same candidates
  cp_keys <- candidates_cp$locus_key
  simple_keys <- candidates_simple$locus_key
  
  common_candidates <- intersect(cp_keys, simple_keys)
  cp_only <- setdiff(cp_keys, simple_keys)
  simple_only <- setdiff(simple_keys, cp_keys)
  
  cat("Common candidates:", length(common_candidates), "\n")
  cat("cp_enrich only:", length(cp_only), "\n")
  cat("simple only:", length(simple_only), "\n")
  
  if (length(cp_only) > 0) {
    cat("cp_enrich only candidates:", paste(head(cp_only, 5), collapse = ", "), "\n")
  }
  if (length(simple_only) > 0) {
    cat("simple only candidates:", paste(head(simple_only, 5), collapse = ", "), "\n")
  }
}

cat("\n")

# Function to analyze gene sets for GO
analyze_go_gene_sets <- function(con, db_name) {
  cat("GO Gene Set Analysis -", db_name, "\n")
  cat(paste(rep("-", 40), collapse = ""), "\n")
  
  # Try to replicate the exact gene set extraction used in ORA
  tryCatch({
    # Use the same function that ORA uses
    go_data <- extract_go_terms_for_enrichment(con, "stored", verbose = TRUE)
    
    cat("  Foreground genes (candidates with GO):", length(go_data$foreground$genes), "\n")
    cat("  Background genes (non-candidates with GO):", length(go_data$background$genes), "\n")
    cat("  Total GO terms available:", nrow(go_data$all_go_terms), "\n")
    
    # Show some example genes
    cat("  Example foreground genes:", paste(head(go_data$foreground$genes, 3), collapse = ", "), "\n")
    cat("  Example background genes:", paste(head(go_data$background$genes, 3), collapse = ", "), "\n")
    
    return(go_data)
    
  }, error = function(e) {
    cat("  ERROR extracting GO data:", e$message, "\n")
    return(NULL)
  })
}

# Function to analyze gene sets for KEGG  
analyze_kegg_gene_sets <- function(con, db_name) {
  cat("KEGG Gene Set Analysis -", db_name, "\n")
  cat(paste(rep("-", 40), collapse = ""), "\n")
  
  tryCatch({
    # Use the same function that ORA uses
    kegg_data <- extract_kegg_terms_for_enrichment(con, "stored", verbose = TRUE)
    
    cat("  Foreground genes (candidates with KEGG):", length(kegg_data$foreground$genes), "\n")
    cat("  Background genes (non-candidates with KEGG):", length(kegg_data$background$genes), "\n")
    cat("  Total KEGG pathways available:", nrow(kegg_data$all_pathways), "\n")
    
    # Show some example genes
    cat("  Example foreground genes:", paste(head(kegg_data$foreground$genes, 3), collapse = ", "), "\n")
    cat("  Example background genes:", paste(head(kegg_data$background$genes, 3), collapse = ", "), "\n")
    
    return(kegg_data)
    
  }, error = function(e) {
    cat("  ERROR extracting KEGG data:", e$message, "\n")
    return(NULL)
  })
}

# Analyze gene sets
cat("4. GENE SET COMPOSITION ANALYSIS\n")
cat("==================================\n")

# GO analysis
go_data_cp <- analyze_go_gene_sets(con_cp, "cp_enrich")
cat("\n")
go_data_simple <- analyze_go_gene_sets(con_simple, "simple")
cat("\n")

# KEGG analysis
kegg_data_cp <- analyze_kegg_gene_sets(con_cp, "cp_enrich")
cat("\n")
kegg_data_simple <- analyze_kegg_gene_sets(con_simple, "simple")
cat("\n")

# Compare gene sets if both successful
if (!is.null(go_data_cp) && !is.null(go_data_simple)) {
  cat("5. GO GENE SET COMPARISON\n")
  cat("===========================\n")
  
  # Compare foreground sets
  fg_cp <- go_data_cp$foreground$genes
  fg_simple <- go_data_simple$foreground$genes
  
  common_fg <- intersect(fg_cp, fg_simple)
  cp_only_fg <- setdiff(fg_cp, fg_simple)
  simple_only_fg <- setdiff(simple_fg, fg_cp)
  
  cat("Foreground gene comparison:\n")
  cat("  Common genes:", length(common_fg), "\n")
  cat("  cp_enrich only:", length(cp_only_fg), "\n")
  cat("  simple only:", length(simple_only_fg), "\n")
  
  # Compare background sets
  bg_cp <- go_data_cp$background$genes
  bg_simple <- go_data_simple$background$genes
  
  common_bg <- intersect(bg_cp, bg_simple)
  cp_only_bg <- setdiff(bg_cp, bg_simple)
  simple_only_bg <- setdiff(simple_bg, bg_cp)
  
  cat("Background gene comparison:\n")
  cat("  Common genes:", length(common_bg), "\n")
  cat("  cp_enrich only:", length(cp_only_bg), "\n")  
  cat("  simple only:", length(simple_only_bg), "\n")
  
  if (length(cp_only_fg) > 0) {
    cat("  cp_enrich only foreground:", paste(head(cp_only_fg, 3), collapse = ", "), "\n")
  }
  if (length(simple_only_fg) > 0) {
    cat("  simple only foreground:", paste(head(simple_only_fg, 3), collapse = ", "), "\n")
  }
  
  cat("\n")
}

# Analyze ORA results if available
cat("6. ORA RESULTS COMPARISON\n")
cat("===========================\n")

# Check ora_analyses table
analyze_ora_results <- function(con, db_name) {
  tables <- DBI::dbListTables(con)
  
  if ("ora_analyses" %in% tables) {
    analyses <- DBI::dbGetQuery(con, "
      SELECT analysis_id, annotation_type, term_type, 
             total_foreground_genes, total_background_genes,
             analysis_date, enrichment_method
      FROM ora_analyses 
      ORDER BY analysis_id
    ")
    
    cat(db_name, "ORA analyses:\n")
    if (nrow(analyses) > 0) {
      print(analyses)
      cat("\n")
    } else {
      cat("  No ORA analyses found\n\n")
    }
    
    return(analyses)
  } else {
    cat(db_name, "- No ora_analyses table\n\n")
    return(data.frame())
  }
}

analyses_cp <- analyze_ora_results(con_cp, "cp_enrich")
analyses_simple <- analyze_ora_results(con_simple, "simple")

# Summary and assessment
cat("7. SUMMARY AND ASSESSMENT\n")
cat("===========================\n")

cat("Key Findings:\n")
cat("-------------\n")

# Check the suspected issue
if (!is.null(go_data_cp) && !is.null(go_data_simple)) {
  cp_bg_size <- length(go_data_cp$background$genes)
  simple_bg_size <- length(go_data_simple$background$genes)
  
  cat("• Background set sizes (GO):\n")
  cat("  - cp_enrich:", cp_bg_size, "genes\n")
  cat("  - simple:", simple_bg_size, "genes\n")
  
  if (abs(cp_bg_size - simple_bg_size) > 10) {
    cat("  ⚠️  SIGNIFICANT DIFFERENCE in background set sizes!\n")
    cat("      This could explain the different enrichment results.\n")
  } else {
    cat("  ✓ Background set sizes are similar\n")
  }
}

cat("• Total annotations:\n")
cat("  - cp_enrich:", stats_cp$total_annotations, "\n")
cat("  - simple:", stats_simple$total_annotations, "\n")

if (stats_cp$total_annotations != stats_simple$total_annotations) {
  cat("  ⚠️  Different total annotation counts!\n")
} else {
  cat("  ✓ Same total annotation counts\n")
}

cat("• Annotations with GO terms:\n")
cat("  - cp_enrich:", stats_cp$go_annotations, "\n") 
cat("  - simple:", stats_simple$go_annotations, "\n")

if (stats_cp$go_annotations != stats_simple$go_annotations) {
  cat("  ⚠️  Different GO annotation counts!\n")
} else {
  cat("  ✓ Same GO annotation counts\n")
}

cat("\nMethodological Assessment:\n")
cat("--------------------------\n")
cat("For ORA background sets, the correct approach is to include:\n")
cat("• ONLY genes that have the relevant annotation type (GO/KEGG)\n")
cat("• This ensures the statistical test is fair and meaningful\n")
cat("• Including genes without annotations dilutes the background\n")
cat("• This can lead to false positives or missed enrichments\n")

cat("\n=== ANALYSIS COMPLETE ===\n")

# Close connections
close_funseq_db(con_cp)
close_funseq_db(con_simple)