# GO Evidence Codes Quick Reference Guide
# 
# This script provides a quick reference for GO evidence codes and recommended
# settings for the annotate_blast_results() function

library(funseqR)

cat("=== GO Evidence Codes Quick Reference ===\n\n")

cat("EXPERIMENTAL EVIDENCE (Highest Confidence):\n")
cat("===========================================\n")
cat("EXP - Inferred from Experiment: Direct experimental evidence\n")
cat("IDA - Inferred from Direct Assay: Direct molecular interaction/localization\n")
cat("IPI - Inferred from Physical Interaction: Protein-protein interactions\n")
cat("IMP - Inferred from Mutant Phenotype: Loss/gain of function studies\n")
cat("IGI - Inferred from Genetic Interaction: Epistasis or suppression\n")
cat("IEP - Inferred from Expression Pattern: Temporal/spatial expression\n\n")

cat("CURATED EVIDENCE (High Confidence):\n")
cat("===================================\n")
cat("TAS - Traceable Author Statement: Curator judgment from literature\n")
cat("IC  - Inferred by Curator: Expert biocuration\n\n")

cat("COMPUTATIONAL EVIDENCE (Medium Confidence):\n")
cat("===========================================\n")
cat("ISS - Inferred from Sequence Similarity: Homology-based transfer\n")
cat("ISO - Inferred from Sequence Orthology: Ortholog-based transfer\n")
cat("ISA - Inferred from Sequence Alignment: Alignment-based transfer\n")
cat("ISM - Inferred from Sequence Model: Domain/motif-based prediction\n")
cat("IGC - Inferred from Genomic Context: Synteny or gene neighborhood\n")
cat("IBA - Inferred from Biological Aspect of Ancestor: Phylogenetic inference\n")
cat("IBD - Inferred from Biological Aspect of Descendant: Phylogenetic inference\n")
cat("IKR - Inferred from Key Residues: Critical amino acid analysis\n")
cat("IRD - Inferred from Rapid Divergence: Evolutionary rate analysis\n")
cat("RCA - Inferred from Reviewed Computational Analysis: Manual review of computational prediction\n\n")

cat("ELECTRONIC EVIDENCE (Lower Confidence):\n")
cat("=======================================\n")
cat("IEA - Inferred from Electronic Annotation: Computational prediction\n\n")

cat("=== RECOMMENDED SETTINGS FOR NON-MODEL ORGANISMS ===\n\n")

cat("1. CONSERVATIVE (High confidence only):\n")
cat("   evidence_keep = c('EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC')\n")
cat("   - Use when you need high-confidence annotations only\n")
cat("   - May result in fewer GO terms\n")
cat("   - Good for functional validation studies\n\n")

cat("2. BALANCED (Recommended for most analyses):\n")
cat("   evidence_keep = c('EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC', 'IEA', 'ISS')\n")
cat("   - Includes computational predictions and homology transfers\n")
cat("   - Good balance of coverage and quality\n")
cat("   - Recommended for GO enrichment analyses\n\n")

cat("3. COMPREHENSIVE (Maximum coverage):\n")
cat("   evidence_keep = NULL\n")
cat("   - Accepts all evidence types\n")
cat("   - Use when functional coverage is more important than evidence quality\n")
cat("   - Good for exploratory analyses\n\n")

cat("=== EXAMPLE USAGE ===\n\n")

cat("# Conservative annotation\n")
cat("annotation_results <- annotate_blast_results(con, blast_param_id,\n")
cat("  evidence_keep = c('EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC'))\n\n")

cat("# Balanced annotation (recommended)\n")
cat("annotation_results <- annotate_blast_results(con, blast_param_id,\n")
cat("  evidence_keep = c('EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC', 'IEA', 'ISS'))\n\n")

cat("# Comprehensive annotation\n")
cat("annotation_results <- annotate_blast_results(con, blast_param_id,\n")
cat("  evidence_keep = NULL)\n\n")

cat("=== COMPOUND EVIDENCE CODES ===\n\n")
cat("UniProt often provides compound evidence codes like:\n")
cat("- IEA:UniProtKB-SubCell (IEA from UniProt subcellular localization)\n")
cat("- IMP:ZFIN (IMP from ZFIN database)\n")
cat("- ISS:UniProtKB (ISS from UniProt)\n\n")
cat("The funseqR package automatically extracts the primary evidence code\n")
cat("(the part before the colon) for filtering purposes.\n\n")

cat("=== CHECKING YOUR RESULTS ===\n\n")
cat("After running annotation, check your results:\n\n")
cat("str(annotation_results)\n")
cat("print(paste('GO terms extracted:', annotation_results$go_terms))\n")
cat("print(paste('Missing GO terms:', annotation_results$verification$missing_go))\n\n")

cat("If you're getting 0 GO terms, try:\n")
cat("1. Include IEA evidence codes\n")
cat("2. Use evidence_keep = NULL for debugging\n")
cat("3. Use debug_accessions to examine specific proteins\n")
cat("4. Check the funseqR documentation: ?annotate_blast_results\n")