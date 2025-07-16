# Investigate ORA differences between cp_enrich and simple databases
# This script will systematically compare the two databases to identify
# differences in gene set composition that could explain ORA result differences

library(DBI)
library(RSQLite)

# Connect to both databases
cp_enrich_db <- dbConnect(SQLite(), "cp_enrich_funseq_project.db")
simple_db <- dbConnect(SQLite(), "simple_funseq_project.db")

cat("=== DATABASE INVESTIGATION: ORA GENE SET COMPOSITION DIFFERENCES ===\n\n")

# 1. DATABASE CONNECTION AND BASIC STATS
cat("1. DATABASE CONNECTION AND BASIC STATS\n")
cat("=====================================\n\n")

# Check tables in both databases
cat("Tables in cp_enrich database:\n")
cp_tables <- dbListTables(cp_enrich_db)
print(cp_tables)

cat("\nTables in simple database:\n")
simple_tables <- dbListTables(simple_db)
print(simple_tables)

# Basic annotation counts
cat("\nBasic Annotation Statistics:\n")
cat("----------------------------\n")

# Total annotations in each database
cp_total_annotations <- dbGetQuery(cp_enrich_db, "SELECT COUNT(*) as count FROM annotations")$count
simple_total_annotations <- dbGetQuery(simple_db, "SELECT COUNT(*) as count FROM annotations")$count

cat("Total annotations:\n")
cat(sprintf("  cp_enrich: %d\n", cp_total_annotations))
cat(sprintf("  simple: %d\n", simple_total_annotations))

# Annotations with GO terms
cp_go_annotations <- dbGetQuery(cp_enrich_db, "SELECT COUNT(DISTINCT locus_id) as count FROM annotations WHERE go_id IS NOT NULL AND go_id != ''")$count
simple_go_annotations <- dbGetQuery(simple_db, "SELECT COUNT(DISTINCT locus_id) as count FROM annotations WHERE go_id IS NOT NULL AND go_id != ''")$count

cat("Unique loci with GO annotations:\n")
cat(sprintf("  cp_enrich: %d\n", cp_go_annotations))
cat(sprintf("  simple: %d\n", simple_go_annotations))

# Annotations with KEGG terms
cp_kegg_annotations <- dbGetQuery(cp_enrich_db, "SELECT COUNT(DISTINCT locus_id) as count FROM annotations WHERE kegg_pathway IS NOT NULL AND kegg_pathway != ''")$count
simple_kegg_annotations <- dbGetQuery(simple_db, "SELECT COUNT(DISTINCT locus_id) as count FROM annotations WHERE kegg_pathway IS NOT NULL AND kegg_pathway != ''")$count

cat("Unique loci with KEGG annotations:\n")
cat(sprintf("  cp_enrich: %d\n", cp_kegg_annotations))
cat(sprintf("  simple: %d\n", simple_kegg_annotations))

cat("\n")

# 2. CANDIDATE LOCI COMPARISON
cat("2. CANDIDATE LOCI COMPARISON\n")
cat("============================\n\n")

# Check if candidate_loci table exists and compare
if ("candidate_loci" %in% cp_tables && "candidate_loci" %in% simple_tables) {
  cp_candidates <- dbGetQuery(cp_enrich_db, "SELECT COUNT(*) as count FROM candidate_loci WHERE is_candidate = 1")$count
  simple_candidates <- dbGetQuery(simple_db, "SELECT COUNT(*) as count FROM candidate_loci WHERE is_candidate = 1")$count
  
  cat("Candidate loci counts:\n")
  cat(sprintf("  cp_enrich: %d\n", cp_candidates))
  cat(sprintf("  simple: %d\n", simple_candidates))
  
  # Check if same loci are candidates
  cp_candidate_ids <- dbGetQuery(cp_enrich_db, "SELECT locus_id FROM candidate_loci WHERE is_candidate = 1")$locus_id
  simple_candidate_ids <- dbGetQuery(simple_db, "SELECT locus_id FROM candidate_loci WHERE is_candidate = 1")$locus_id
  
  shared_candidates <- length(intersect(cp_candidate_ids, simple_candidate_ids))
  cp_only <- length(setdiff(cp_candidate_ids, simple_candidate_ids))
  simple_only <- length(setdiff(simple_candidate_ids, cp_candidate_ids))
  
  cat("Candidate loci overlap:\n")
  cat(sprintf("  Shared candidates: %d\n", shared_candidates))
  cat(sprintf("  cp_enrich only: %d\n", cp_only))
  cat(sprintf("  simple only: %d\n", simple_only))
} else {
  cat("Candidate loci table not found in one or both databases\n")
}

cat("\n")

# 3. GENE SET ANALYSIS - GO
cat("3. GENE SET ANALYSIS - GO\n")
cat("=========================\n\n")

# Check for stored_candidate_loci table (newer approach)
if ("stored_candidate_loci" %in% cp_tables) {
  cat("Using stored_candidate_loci table for cp_enrich database\n")
  
  # Get GO gene sets for cp_enrich
  cp_go_foreground <- dbGetQuery(cp_enrich_db, "
    SELECT DISTINCT a.locus_id 
    FROM annotations a 
    JOIN stored_candidate_loci scl ON a.locus_id = scl.locus_id 
    WHERE scl.is_candidate = 1 
    AND a.go_id IS NOT NULL 
    AND a.go_id != ''
  ")$locus_id
  
  cp_go_background <- dbGetQuery(cp_enrich_db, "
    SELECT DISTINCT a.locus_id 
    FROM annotations a 
    JOIN stored_candidate_loci scl ON a.locus_id = scl.locus_id 
    WHERE scl.is_candidate = 0 
    AND a.go_id IS NOT NULL 
    AND a.go_id != ''
  ")$locus_id
  
} else if ("candidate_loci" %in% cp_tables) {
  cat("Using candidate_loci table for cp_enrich database\n")
  
  # Get GO gene sets for cp_enrich
  cp_go_foreground <- dbGetQuery(cp_enrich_db, "
    SELECT DISTINCT a.locus_id 
    FROM annotations a 
    JOIN candidate_loci cl ON a.locus_id = cl.locus_id 
    WHERE cl.is_candidate = 1 
    AND a.go_id IS NOT NULL 
    AND a.go_id != ''
  ")$locus_id
  
  cp_go_background <- dbGetQuery(cp_enrich_db, "
    SELECT DISTINCT a.locus_id 
    FROM annotations a 
    JOIN candidate_loci cl ON a.locus_id = cl.locus_id 
    WHERE cl.is_candidate = 0 
    AND a.go_id IS NOT NULL 
    AND a.go_id != ''
  ")$locus_id
  
} else {
  cat("No candidate loci table found in cp_enrich database\n")
  cp_go_foreground <- c()
  cp_go_background <- c()
}

# Get GO gene sets for simple database
if ("candidate_loci" %in% simple_tables) {
  cat("Using candidate_loci table for simple database\n")
  
  simple_go_foreground <- dbGetQuery(simple_db, "
    SELECT DISTINCT a.locus_id 
    FROM annotations a 
    JOIN candidate_loci cl ON a.locus_id = cl.locus_id 
    WHERE cl.is_candidate = 1 
    AND a.go_id IS NOT NULL 
    AND a.go_id != ''
  ")$locus_id
  
  simple_go_background <- dbGetQuery(simple_db, "
    SELECT DISTINCT a.locus_id 
    FROM annotations a 
    JOIN candidate_loci cl ON a.locus_id = cl.locus_id 
    WHERE cl.is_candidate = 0 
    AND a.go_id IS NOT NULL 
    AND a.go_id != ''
  ")$locus_id
  
} else {
  cat("No candidate loci table found in simple database\n")
  simple_go_foreground <- c()
  simple_go_background <- c()
}

cat("GO Gene Set Sizes:\n")
cat(sprintf("  cp_enrich foreground: %d\n", length(cp_go_foreground)))
cat(sprintf("  cp_enrich background: %d\n", length(cp_go_background)))
cat(sprintf("  simple foreground: %d\n", length(simple_go_foreground)))
cat(sprintf("  simple background: %d\n", length(simple_go_background)))

# Check if background includes genes without GO terms
cp_total_background <- dbGetQuery(cp_enrich_db, "
  SELECT COUNT(DISTINCT locus_id) as count 
  FROM candidate_loci 
  WHERE is_candidate = 0
")$count

simple_total_background <- dbGetQuery(simple_db, "
  SELECT COUNT(DISTINCT locus_id) as count 
  FROM candidate_loci 
  WHERE is_candidate = 0
")$count

cat("Total background loci (including those without GO):\n")
cat(sprintf("  cp_enrich: %d\n", cp_total_background))
cat(sprintf("  simple: %d\n", simple_total_background))

cat("\n")

# 4. GENE SET ANALYSIS - KEGG
cat("4. GENE SET ANALYSIS - KEGG\n")
cat("===========================\n\n")

# Similar analysis for KEGG
if ("stored_candidate_loci" %in% cp_tables) {
  cp_kegg_foreground <- dbGetQuery(cp_enrich_db, "
    SELECT DISTINCT a.locus_id 
    FROM annotations a 
    JOIN stored_candidate_loci scl ON a.locus_id = scl.locus_id 
    WHERE scl.is_candidate = 1 
    AND a.kegg_pathway IS NOT NULL 
    AND a.kegg_pathway != ''
  ")$locus_id
  
  cp_kegg_background <- dbGetQuery(cp_enrich_db, "
    SELECT DISTINCT a.locus_id 
    FROM annotations a 
    JOIN stored_candidate_loci scl ON a.locus_id = scl.locus_id 
    WHERE scl.is_candidate = 0 
    AND a.kegg_pathway IS NOT NULL 
    AND a.kegg_pathway != ''
  ")$locus_id
  
} else if ("candidate_loci" %in% cp_tables) {
  cp_kegg_foreground <- dbGetQuery(cp_enrich_db, "
    SELECT DISTINCT a.locus_id 
    FROM annotations a 
    JOIN candidate_loci cl ON a.locus_id = cl.locus_id 
    WHERE cl.is_candidate = 1 
    AND a.kegg_pathway IS NOT NULL 
    AND a.kegg_pathway != ''
  ")$locus_id
  
  cp_kegg_background <- dbGetQuery(cp_enrich_db, "
    SELECT DISTINCT a.locus_id 
    FROM annotations a 
    JOIN candidate_loci cl ON a.locus_id = cl.locus_id 
    WHERE cl.is_candidate = 0 
    AND a.kegg_pathway IS NOT NULL 
    AND a.kegg_pathway != ''
  ")$locus_id
  
} else {
  cp_kegg_foreground <- c()
  cp_kegg_background <- c()
}

if ("candidate_loci" %in% simple_tables) {
  simple_kegg_foreground <- dbGetQuery(simple_db, "
    SELECT DISTINCT a.locus_id 
    FROM annotations a 
    JOIN candidate_loci cl ON a.locus_id = cl.locus_id 
    WHERE cl.is_candidate = 1 
    AND a.kegg_pathway IS NOT NULL 
    AND a.kegg_pathway != ''
  ")$locus_id
  
  simple_kegg_background <- dbGetQuery(simple_db, "
    SELECT DISTINCT a.locus_id 
    FROM annotations a 
    JOIN candidate_loci cl ON a.locus_id = cl.locus_id 
    WHERE cl.is_candidate = 0 
    AND a.kegg_pathway IS NOT NULL 
    AND a.kegg_pathway != ''
  ")$locus_id
  
} else {
  simple_kegg_foreground <- c()
  simple_kegg_background <- c()
}

cat("KEGG Gene Set Sizes:\n")
cat(sprintf("  cp_enrich foreground: %d\n", length(cp_kegg_foreground)))
cat(sprintf("  cp_enrich background: %d\n", length(cp_kegg_background)))
cat(sprintf("  simple foreground: %d\n", length(simple_kegg_foreground)))
cat(sprintf("  simple background: %d\n", length(simple_kegg_background)))

cat("\n")

# 5. ORA RESULTS COMPARISON
cat("5. ORA RESULTS COMPARISON\n")
cat("=========================\n\n")

# Check ora_analyses table
if ("ora_analyses" %in% cp_tables && "ora_analyses" %in% simple_tables) {
  cat("ORA Analyses parameters:\n")
  
  cp_ora <- dbGetQuery(cp_enrich_db, "SELECT * FROM ora_analyses ORDER BY id DESC LIMIT 5")
  simple_ora <- dbGetQuery(simple_db, "SELECT * FROM ora_analyses ORDER BY id DESC LIMIT 5")
  
  cat("cp_enrich recent analyses:\n")
  print(cp_ora)
  
  cat("\nsimple recent analyses:\n")
  print(simple_ora)
  
  # Check ora_results
  if ("ora_results" %in% cp_tables && "ora_results" %in% simple_tables) {
    cp_results_count <- dbGetQuery(cp_enrich_db, "SELECT COUNT(*) as count FROM ora_results")$count
    simple_results_count <- dbGetQuery(simple_db, "SELECT COUNT(*) as count FROM ora_results")$count
    
    cat(sprintf("\nORA results count:\n"))
    cat(sprintf("  cp_enrich: %d\n", cp_results_count))
    cat(sprintf("  simple: %d\n", simple_results_count))
    
    # Sample results
    if (cp_results_count > 0) {
      cat("\ncp_enrich sample results (most significant):\n")
      cp_sample <- dbGetQuery(cp_enrich_db, "SELECT * FROM ora_results ORDER BY p_value LIMIT 5")
      print(cp_sample)
    }
    
    if (simple_results_count > 0) {
      cat("\nsimple sample results (most significant):\n")
      simple_sample <- dbGetQuery(simple_db, "SELECT * FROM ora_results ORDER BY p_value LIMIT 5")
      print(simple_sample)
    }
  }
}

cat("\n")

# 6. DETAILED BACKGROUND INVESTIGATION
cat("6. DETAILED BACKGROUND INVESTIGATION\n")
cat("====================================\n\n")

# Check if background in either database includes ALL loci or only annotated loci
cat("Background composition analysis:\n")

# For cp_enrich database
if ("stored_candidate_loci" %in% cp_tables || "candidate_loci" %in% cp_tables) {
  # All loci in the candidate table
  cp_all_loci <- dbGetQuery(cp_enrich_db, "SELECT COUNT(DISTINCT locus_id) as count FROM candidate_loci")$count
  
  # Background loci with any annotation
  cp_bg_with_annotations <- dbGetQuery(cp_enrich_db, "
    SELECT COUNT(DISTINCT cl.locus_id) as count 
    FROM candidate_loci cl 
    JOIN annotations a ON cl.locus_id = a.locus_id 
    WHERE cl.is_candidate = 0
  ")$count
  
  # Background loci without any annotations
  cp_bg_without_annotations <- cp_total_background - cp_bg_with_annotations
  
  cat("cp_enrich background composition:\n")
  cat(sprintf("  Total background loci: %d\n", cp_total_background))
  cat(sprintf("  Background with annotations: %d\n", cp_bg_with_annotations))
  cat(sprintf("  Background without annotations: %d\n", cp_bg_without_annotations))
  
  # Check if unannotated loci are being included in GO analysis
  cp_bg_no_go <- cp_total_background - length(cp_go_background)
  cat(sprintf("  Background loci without GO terms: %d\n", cp_bg_no_go))
}

# For simple database
if ("candidate_loci" %in% simple_tables) {
  # Background loci with any annotation
  simple_bg_with_annotations <- dbGetQuery(simple_db, "
    SELECT COUNT(DISTINCT cl.locus_id) as count 
    FROM candidate_loci cl 
    JOIN annotations a ON cl.locus_id = a.locus_id 
    WHERE cl.is_candidate = 0
  ")$count
  
  # Background loci without any annotations
  simple_bg_without_annotations <- simple_total_background - simple_bg_with_annotations
  
  cat("simple background composition:\n")
  cat(sprintf("  Total background loci: %d\n", simple_total_background))
  cat(sprintf("  Background with annotations: %d\n", simple_bg_with_annotations))
  cat(sprintf("  Background without annotations: %d\n", simple_bg_without_annotations))
  
  # Check if unannotated loci are being included in GO analysis
  simple_bg_no_go <- simple_total_background - length(simple_go_background)
  cat(sprintf("  Background loci without GO terms: %d\n", simple_bg_no_go))
}

cat("\n")

# Close database connections
dbDisconnect(cp_enrich_db)
dbDisconnect(simple_db)

cat("=== INVESTIGATION COMPLETE ===\n")