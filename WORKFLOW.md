# funseqR Core Workflow Documentation

## Overview

### Purpose and Scope

The funseqR package addresses a critical gap in functional genomics tools for **non-model species** research. While numerous packages exist for functional annotation and enrichment analysis, most assume access to well-curated genome annotations, established gene models, and comprehensive functional databases that are typically only available for major model organisms (human, mouse, *Drosophila*, *Arabidopsis*, etc.).

### Target Applications

funseqR is specifically designed for **population and conservation genomics studies** that work with:

- **Large-scale genomic variant datasets**: 10,000s to 100,000s of SNPs
- **Reduced-representation sequencing**: ddRAD-seq, RAD-seq, GBS data
- **Whole genome resequencing**: Population-scale variant discovery
- **Candidate adaptive loci**: SNPs identified through genotype-environment association (GEA) analyses
- **Outlier loci**: Variants showing signatures of selection or local adaptation
- **Conservation genomics**: Functional characterization of adaptive variation

### Filling the Non-Model Species Gap

**Traditional functional genomics packages** (e.g., DAVID, Enrichr, clusterProfiler, topGO) are optimized for:
- Model organisms with complete genome annotations
- Gene expression data from RNA-seq experiments  
- Predefined gene sets and pathway databases
- Well-characterized protein-coding genes with established functional annotations

**funseqR addresses the challenge** when working with:
- **Non-model species** lacking comprehensive genome annotations
- **Anonymous genomic loci** without known gene associations
- **Population-scale variant data** rather than expression data
- **Conservation genomics contexts** where functional insights are critical for management decisions

### Methodological Innovation

Rather than relying on pre-existing annotations, funseqR takes a **discovery-based approach**:

1. **Direct sequence analysis**: Uses BLAST to functionally characterize genomic regions surrounding variants
2. **Protein database integration**: Leverages the comprehensive UniProt database rather than species-specific annotations
3. **Flexible variant input**: Works with any genomic variant regardless of annotation status
4. **Comparative enrichment**: Enables functional comparison between candidate and background loci sets

### Research Context

This approach is particularly valuable for studies investigating:

- **Local adaptation**: Functional characterization of environmentally-associated genetic variants
- **Population differentiation**: Understanding the functional basis of population structure
- **Conservation genomics**: Identifying functionally important regions for conservation priorities
- **Evolutionary genomics**: Characterizing the genetic basis of adaptive traits
- **Climate change research**: Functional annotation of climate-associated genetic variants

### Contrast with Existing Tools

| **Traditional Tools** | **funseqR** |
|---|---|
| Model organism focus | Non-model species emphasis |
| Gene expression input | Genomic variant input |
| Predefined gene sets | Discovery-based annotation |
| Species-specific databases | Universal protein databases |
| Well-annotated genomes | Anonymous genomic loci |
| Laboratory/clinical studies | Field-based ecological studies |

The funseqR package provides a complete workflow for functional annotation of genomic variants through BLAST searching and UniProt annotation, specifically designed for the unique challenges of non-model species research. This document outlines each step in the core workflow, the computational approach used, and the outputs generated.

## Workflow Architecture

The funseqR workflow follows a **one-project-per-database** design where each SQLite database represents a single analysis project. Data flows through the following stages:

```
VCF Files → Database Import → Flanking Sequence Extraction → BLAST Search → Functional Annotation → GO Enrichment
```

---

## Step 1: Database Creation

### Function: `create_funseq_db()`

**Purpose**: Initialize a new SQLite database with the complete funseqR schema.

**Approach**:
- Creates 12 core tables with proper foreign key relationships
- Establishes indexes for efficient querying
- Includes tables for: input files, VCF data, reference genomes, sequences, BLAST parameters/results, annotations, GO terms, KEGG references, and UniProt cache

**Output**:
```r
con <- create_funseq_db("analysis.db", force = TRUE)
# Returns: Database connection object
```

**Database Schema Created**:
- `input_files` - File registry and metadata
- `vcf_data` - Variant call data  
- `reference_genomes` - Genome metadata
- `reference_sequences` - Sequence data (optional storage)
- `flanking_sequences` - Extracted sequences around variants
- `blast_parameters` - BLAST search configurations
- `blast_results` - BLAST hit results
- `annotations` - UniProt protein annotations
- `go_terms` - Gene Ontology associations
- `kegg_references` - KEGG pathway references
- `uniprot_cache` - API response caching
- `go_enrichment_analyses` / `go_enrichment_results` - Enrichment analysis storage

---

## Step 2: Data Import

### VCF Import: `import_vcf_to_db()`

**Purpose**: Import variant call format (VCF) files into the database.

**Approach**:
- Registers file metadata (path, hash, import date)
- Parses VCF using `vcfR` package
- Extracts core variant information: chromosome, position, reference/alternative alleles
- Stores variant-level metadata (quality scores, filters, INFO fields)

**Input**: VCF file path
**Output**:
```r
vcf_import <- import_vcf_to_db(con, "variants.vcf")
# Returns: list(file_id = 1, vcf_count = 14289)
```

**Data Stored**:
- File registry entry with SHA256 hash
- Individual variant records with genomic coordinates
- Quality metrics and variant annotations

### Reference Genome Import: `import_reference_to_db()`

**Purpose**: Import reference genome sequences for flanking sequence extraction.

**Approach**:
- Registers FASTA file metadata
- Reads sequences using `Biostrings::readDNAStringSet()`
- Optionally stores sequences in database or maintains file references
- Creates sequence index for rapid coordinate-based retrieval

**Input**: FASTA file path, genome name, build version
**Output**:
```r
ref_import <- import_reference_to_db(con, "genome.fasta", "Species_name", "v1.0")
# Returns: list(file_id = 2, genome_id = 1, sequence_count = 15234)
```

**Data Stored**:
- Genome metadata (name, build, source file)
- Sequence records (name, length, optional sequence data)
- Coordinate indexing for efficient retrieval

---

## Step 3: Flanking Sequence Extraction

### Function: `import_flanking_seqs_to_db()`

**Purpose**: Extract DNA sequences surrounding each variant for BLAST searching.

**Approach**:
- For each variant, extracts configurable flanking regions (default: ±500bp)
- Optional ORF extraction for more targeted searches
- Handles chromosome boundaries and sequence availability
- Stores sequences with coordinate metadata for traceability
- Uses `Biostrings::subseq()` for efficient sequence extraction

**Input**: VCF file ID, genome ID, flank size, optional ORF parameters
**Output**:
```r
# Raw sequences only
flanking_import <- import_flanking_seqs_to_db(con, vcf_import$file_id, ref_import$genome_id, flank_size = 500)
# Returns: list(flanking_count = 14289, orf_count = 0)

# With ORF extraction
flanking_import <- import_flanking_seqs_to_db(con, vcf_import$file_id, ref_import$genome_id, 
                                              flank_size = 500, translate_flanks = TRUE)
# Returns: list(flanking_count = 14289, orf_count = 8234)
```

**Data Stored**:
- Flanking sequence records linked to variants
- Coordinate information (start/end positions)
- Actual DNA sequences for BLAST input

**Quality Considerations**:
- Sequences containing 'N' characters are flagged
- Minimum sequence length requirements
- Chromosome boundary handling

---

## Step 4: BLAST Search

### Function: `perform_blast_db()`

**Purpose**: Search flanking sequences against protein databases to identify potential coding regions.

**Approach**:
- Registers search parameters for reproducibility
- Exports flanking sequences to temporary FASTA files
- Executes BLAST or DIAMOND via system calls with configurable parameters
- Parses tabular output (-outfmt 6)
- Stores results with full provenance tracking

**Input**: VCF file ID, database path/name, search parameters, engine type
**Output**:
```r
# Traditional BLAST
blast_results <- perform_blast_db(con, vcf_import$file_id, "/path/to/db", "protein_db", "blastx")
# Returns: list(blast_param_id = 1, result_count = 31777, output_base = "...", metadata_id = 1)

# DIAMOND (20-50x faster)
blast_results <- perform_blast_db(con, vcf_import$file_id, "/path/to/db", "protein_db", "diamond_blastx")
# Returns: list(blast_param_id = 2, result_count = 31777, output_base = "...", metadata_id = 2)
```

**Search Configuration**:
- **Types**: blastn/blastx (NCBI BLAST) or diamond_blastn/diamond_blastx (DIAMOND)
- **E-value threshold**: 1e-5 (configurable)
- **Max hits per query**: 5 (configurable)
- **Output format**: Tabular with standard fields
- **Threading**: Configurable CPU usage
- **Sequence types**: raw, orf_nuc, or orf_aa

**Data Stored**:
- BLAST parameter records (reproducibility)
- Individual hit results with alignment statistics
- Database metadata (version, composition, date)
- Hit descriptions and accession numbers

**Quality Metrics**:
- E-value significance thresholds
- Bit score rankings
- Percent identity cutoffs
- Alignment length considerations

---

## Step 5: Functional Annotation

### Function: `annotate_blast_results()`

**Purpose**: Retrieve detailed functional annotations for BLAST hits from UniProt.

**Approach**:
- Extracts unique protein accessions from BLAST results
- Implements multi-tier caching system:
  1. Database cache (fastest)
  2. Live API queries (with rate limiting)
- Parses UniProt JSON responses for functional data
- Stores structured annotations with full provenance

**Input**: BLAST parameter ID, annotation thresholds
**Output**:
```r
annotation_results <- annotate_blast_results(con, blast_results$blast_param_id, max_hits = 5, e_value_threshold = 1e-10)
# Returns: list(annotation_count = 17269, cache_hits = 1205, api_calls = 3847)
```

**Annotation Data Retrieved**:
- **Protein Information**: Accession, entry name, gene names
- **Gene Ontology**: GO IDs, terms, categories (BP/MF/CC), evidence codes
- **KEGG Pathways**: Pathway IDs and names
- **Functional Descriptions**: Protein names and functions

**Caching Strategy**:
- Persistent SQLite cache for offline analysis
- Rate-limited API calls (1-2 seconds between requests)
- JSON response storage for data provenance
- Configurable batch processing

**Quality Control**:
- E-value filtering for high-confidence hits
- UniProt API timeout handling
- Data validation and error recovery
- Annotation completeness tracking

---

## Step 6: Annotation Retrieval

### Function: `get_annotations()`

**Purpose**: Retrieve and consolidate all functional annotations for downstream analysis.

**Approach**:
- Joins BLAST results with protein annotations
- Optionally includes GO terms and KEGG pathway data
- Links annotations back to original genomic variants
- Provides structured output for analysis workflows

**Input**: BLAST parameter ID, data inclusion options
**Output**:
```r
annotations <- get_annotations(con, blast_results$blast_param_id, include_go = TRUE, include_kegg = TRUE)
# Returns: list(annotations = data.frame(17269 obs), go_terms = data.frame(180549 obs), kegg_refs = data.frame(11779 obs))
```

**Consolidated Data Structure**:
- **Primary annotations**: Protein details linked to genomic variants
- **GO terms**: Complete ontology associations with evidence
- **KEGG references**: Pathway associations and functional categories
- **Genomic context**: Original variant coordinates and alleles

---

## Step 7: GO Enrichment Analysis

### Function: `run_go_enrichment_workflow()`

**Purpose**: Perform Gene Ontology enrichment analysis comparing candidate vs background variant sets.

**Approach**:
- Imports candidate variant set and links to background annotations
- Extracts GO term associations for both datasets
- Performs hypergeometric enrichment testing
- Applies multiple testing correction (FDR)
- Generates visualization outputs

**Input**: Candidate VCF file, background file ID, ontology selection
**Output**:
```r
go_results <- run_go_enrichment_workflow(con, "candidates.vcf", background_file_id = vcf_import$file_id)
# Returns: list(status = "success", summary = list(...), enrichment_results = list(...), plots = list(...))
```

**Statistical Analysis**:
- **Test**: Hypergeometric distribution
- **Correction**: False Discovery Rate (FDR)
- **Filtering**: Minimum/maximum gene set sizes
- **Ontologies**: Biological Process, Molecular Function, Cellular Component

**Visualization Outputs**:
- Bubble plots of enriched terms
- Treemap visualizations (if available)
- Comparison plots across ontologies
- Summary tables with statistical metrics

---

## Data Flow Summary

1. **Input**: VCF files + Reference genome
2. **Processing**: Flanking extraction → BLAST search → Annotation
3. **Analysis**: GO enrichment + Statistical testing
4. **Output**: Functional annotations + Enriched pathways + Visualizations

## Database Design Benefits

- **Reproducibility**: Complete parameter and provenance tracking
- **Efficiency**: Intelligent caching and indexing
- **Scalability**: Handles large variant sets efficiently
- **Integration**: Links genomic variants to functional annotations
- **Persistence**: All data permanently stored for reanalysis

## Quality Assurance

- SHA256 file hashing for data integrity
- Parameter tracking for reproducible analysis
- Multi-level validation and error handling
- Configurable quality thresholds throughout pipeline
- Comprehensive logging and progress reporting

---

*This workflow provides a complete pipeline from raw genomic variants to functional insights, with emphasis on reproducibility, data provenance, and statistical rigor.*