-- SQL queries to investigate ORA differences between databases

-- First, examine cp_enrich database
.open cp_enrich_funseq_project.db

.echo on
.headers on

-- Show all tables
.print "=== CP_ENRICH DATABASE TABLES ==="
.tables

-- Basic annotation statistics
.print "\n=== CP_ENRICH BASIC STATISTICS ==="
.print "Total annotations:"
SELECT COUNT(*) as total_annotations FROM annotations;

.print "Unique loci with GO annotations:"
SELECT COUNT(DISTINCT locus_id) as loci_with_go FROM annotations WHERE go_id IS NOT NULL AND go_id != '';

.print "Unique loci with KEGG annotations:"
SELECT COUNT(DISTINCT locus_id) as loci_with_kegg FROM annotations WHERE kegg_pathway IS NOT NULL AND kegg_pathway != '';

-- Check candidate loci
.print "Candidate loci counts:"
SELECT 
  is_candidate,
  COUNT(*) as count
FROM candidate_loci 
GROUP BY is_candidate;

-- Check stored candidate loci if exists
.print "Stored candidate loci (if exists):"
SELECT name FROM sqlite_master WHERE type='table' AND name='stored_candidate_loci';

-- GO foreground/background for candidate_loci
.print "\nGO gene sets using candidate_loci:"
.print "GO foreground (candidates with GO):"
SELECT COUNT(DISTINCT a.locus_id) as go_foreground
FROM annotations a 
JOIN candidate_loci cl ON a.locus_id = cl.locus_id 
WHERE cl.is_candidate = 1 
AND a.go_id IS NOT NULL 
AND a.go_id != '';

.print "GO background (non-candidates with GO):"
SELECT COUNT(DISTINCT a.locus_id) as go_background
FROM annotations a 
JOIN candidate_loci cl ON a.locus_id = cl.locus_id 
WHERE cl.is_candidate = 0 
AND a.go_id IS NOT NULL 
AND a.go_id != '';

-- KEGG foreground/background
.print "\nKEGG gene sets using candidate_loci:"
.print "KEGG foreground (candidates with KEGG):"
SELECT COUNT(DISTINCT a.locus_id) as kegg_foreground
FROM annotations a 
JOIN candidate_loci cl ON a.locus_id = cl.locus_id 
WHERE cl.is_candidate = 1 
AND a.kegg_pathway IS NOT NULL 
AND a.kegg_pathway != '';

.print "KEGG background (non-candidates with KEGG):"
SELECT COUNT(DISTINCT a.locus_id) as kegg_background
FROM annotations a 
JOIN candidate_loci cl ON a.locus_id = cl.locus_id 
WHERE cl.is_candidate = 0 
AND a.kegg_pathway IS NOT NULL 
AND a.kegg_pathway != '';

-- ORA analyses information
.print "\nORA analyses (recent 3):"
SELECT * FROM ora_analyses ORDER BY id DESC LIMIT 3;

.print "\nORA results count:"
SELECT COUNT(*) as total_ora_results FROM ora_results;

.print "\nTop 5 most significant GO results:"
SELECT term_id, term_description, p_value, q_value FROM ora_results 
WHERE analysis_type = 'GO' 
ORDER BY p_value LIMIT 5;

-- Now examine simple database
.close
.open simple_funseq_project.db

.print "\n\n=== SIMPLE DATABASE TABLES ==="
.tables

-- Basic annotation statistics
.print "\n=== SIMPLE BASIC STATISTICS ==="
.print "Total annotations:"
SELECT COUNT(*) as total_annotations FROM annotations;

.print "Unique loci with GO annotations:"
SELECT COUNT(DISTINCT locus_id) as loci_with_go FROM annotations WHERE go_id IS NOT NULL AND go_id != '';

.print "Unique loci with KEGG annotations:"
SELECT COUNT(DISTINCT locus_id) as loci_with_kegg FROM annotations WHERE kegg_pathway IS NOT NULL AND kegg_pathway != '';

-- Check candidate loci
.print "Candidate loci counts:"
SELECT 
  is_candidate,
  COUNT(*) as count
FROM candidate_loci 
GROUP BY is_candidate;

-- GO foreground/background
.print "\nGO gene sets:"
.print "GO foreground (candidates with GO):"
SELECT COUNT(DISTINCT a.locus_id) as go_foreground
FROM annotations a 
JOIN candidate_loci cl ON a.locus_id = cl.locus_id 
WHERE cl.is_candidate = 1 
AND a.go_id IS NOT NULL 
AND a.go_id != '';

.print "GO background (non-candidates with GO):"
SELECT COUNT(DISTINCT a.locus_id) as go_background
FROM annotations a 
JOIN candidate_loci cl ON a.locus_id = cl.locus_id 
WHERE cl.is_candidate = 0 
AND a.go_id IS NOT NULL 
AND a.go_id != '';

-- KEGG foreground/background
.print "\nKEGG gene sets:"
.print "KEGG foreground (candidates with KEGG):"
SELECT COUNT(DISTINCT a.locus_id) as kegg_foreground
FROM annotations a 
JOIN candidate_loci cl ON a.locus_id = cl.locus_id 
WHERE cl.is_candidate = 1 
AND a.kegg_pathway IS NOT NULL 
AND a.kegg_pathway != '';

.print "KEGG background (non-candidates with KEGG):"
SELECT COUNT(DISTINCT a.locus_id) as kegg_background
FROM annotations a 
JOIN candidate_loci cl ON a.locus_id = cl.locus_id 
WHERE cl.is_candidate = 0 
AND a.kegg_pathway IS NOT NULL 
AND a.kegg_pathway != '';

-- ORA analyses information
.print "\nORA analyses (recent 3):"
SELECT * FROM ora_analyses ORDER BY id DESC LIMIT 3;

.print "\nORA results count:"
SELECT COUNT(*) as total_ora_results FROM ora_results;

.print "\nTop 5 most significant GO results:"
SELECT term_id, term_description, p_value, q_value FROM ora_results 
WHERE analysis_type = 'GO' 
ORDER BY p_value LIMIT 5;

.close