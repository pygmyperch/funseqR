# FunseqR Database Schema Documentation

Complete documentation of all 21 tables in the funseqR SQLite database schema.

## Table Summary

| # | Table Name | Purpose | Status | Key Functions |
|---|------------|---------|--------|---------------|
| 1 | metadata | Database versioning & validation | ✅ Keep | `create_funseq_db()`, `connect_funseq_db()` |
| 2 | input_files | File registry & tracking | ⚠️ Legacy | `register_input_file()`, `import_vcf()` |
| 3 | vcf_data | Variant call data storage | ⚠️ Legacy | `import_vcf()` |
| 4 | reference_genomes | Reference genome metadata | ⚠️ Legacy | `import_reference_genome()` |
| 5 | reference_sequences | Reference sequence storage | ⚠️ Legacy | `import_reference_genome()` |
| 6 | flanking_sequences | Extracted sequences around variants | ✅ Keep | `extract_flanking_sequences()` |
| 7 | blast_parameters | BLAST run parameters | 🚫 **REMOVE** | `run_blast()` |
| 8 | blast_results | BLAST search results | ✅ Keep | `run_blast()` |
| 9 | blast_database_metadata | BLAST database info | 🚫 **REMOVE** | `run_blast()` |
| 10 | annotations | UniProt protein annotations | ✅ Keep | `retrieve_annotations()` |
| 11 | go_terms | Gene Ontology terms | ✅ Keep | `retrieve_annotations()` |
| 12 | kegg_references | KEGG pathway references | ✅ Keep | `retrieve_annotations()` |
| 13 | pfam_domains | Pfam domain annotations | ✅ Keep | `retrieve_annotations()` |
| 14 | interpro_families | InterPro family annotations | ✅ Keep | `retrieve_annotations()` |
| 15 | eggnog_categories | eggNOG category annotations | ✅ Keep | `retrieve_annotations()` |
| 16 | analysis_reports | Generated report metadata | ⚠️ Review | Report functions |
| 17 | ora_analyses | Enrichment analysis metadata | ✅ Keep* | `ora()`, `store_ora_results()` |
| 18 | ora_results | Enrichment analysis results | ✅ Keep | `ora()`, `store_ora_results()` |
| 19 | uniprot_cache | UniProt API response cache | ✅ Keep | `retrieve_annotations()` |
| 20 | locus_statistics | Statistical values per locus | ✅ Keep | `define_locus_statistics()` |
| 21 | candidate_loci | Defined candidate loci | ✅ Keep | `define_locus_statistics()` |

*\* = Contains legacy `blast_param_id` column to remove*

---

## Detailed Table Documentation

### 1. metadata
**Purpose:** Database versioning, validation, and creation tracking.

**Schema:**
```sql
CREATE TABLE metadata (
  key TEXT PRIMARY KEY,
  value TEXT NOT NULL
)
```

**What gets stored:**
- `funseqR_version` - Package version used to create database
- `r_version` - R version used to create database  
- `creation_date` - Database creation timestamp

**When populated:** `create_funseq_db()` during database initialization

**How used:** `connect_funseq_db()` validates database and displays version info

**Status:** ✅ **Keep** - Essential for database integrity

---

### 2. input_files ⚠️
**Purpose:** Registry and tracking of imported files (VCF, FASTA, etc.)

**Schema:**
```sql
CREATE TABLE input_files (
  file_id INTEGER PRIMARY KEY,
  file_type TEXT NOT NULL,           -- "vcf", "fasta", etc.
  file_name TEXT NOT NULL,
  file_path TEXT NOT NULL,
  file_hash TEXT,                    -- SHA256 hash for integrity
  import_date TEXT NOT NULL
)
```

**When populated:** 
- `register_input_file()` when any file is imported
- `import_vcf()`, `import_reference_genome()` call register_input_file()

**How used:** Referenced by `vcf_data`, `reference_genomes` via foreign keys

**Status:** ⚠️ **Legacy** - In "one database = one analysis" model, may not need file tracking

---

### 3. vcf_data ⚠️
**Purpose:** Storage of variant call format (VCF) data

**Schema:**
```sql
CREATE TABLE vcf_data (
  vcf_id INTEGER PRIMARY KEY,
  file_id INTEGER NOT NULL,          -- → input_files(file_id)
  chromosome TEXT NOT NULL,
  position INTEGER NOT NULL,
  id TEXT,
  ref TEXT NOT NULL,                 -- Reference allele
  alt TEXT NOT NULL,                 -- Alternative allele
  qual REAL,                         -- Quality score
  filter TEXT,
  info TEXT,
  format TEXT,
  sample_data TEXT,
  FOREIGN KEY (file_id) REFERENCES input_files (file_id)
)
```

**When populated:** `import_vcf()` reads VCF files and stores variants

**How used:** 
- `extract_flanking_sequences()` uses coordinates
- `define_locus_statistics()` can match VCF coordinates to statistics

**Status:** ⚠️ **Legacy** - In streamlined workflow, users provide coordinates directly via `locus_statistics`

---

### 4. reference_genomes ⚠️
**Purpose:** Metadata about imported reference genomes

**Schema:**
```sql
CREATE TABLE reference_genomes (
  genome_id INTEGER PRIMARY KEY,
  file_id INTEGER NOT NULL,          -- → input_files(file_id)
  genome_name TEXT NOT NULL,
  genome_build TEXT,
  FOREIGN KEY (file_id) REFERENCES input_files (file_id)
)
```

**When populated:** `import_reference_genome()` when FASTA files imported

**How used:** Referenced by `reference_sequences` table

**Status:** ⚠️ **Legacy** - May be simplified or removed in streamlined workflow

---

### 5. reference_sequences ⚠️
**Purpose:** Storage of reference genome sequence data

**Schema:**
```sql
CREATE TABLE reference_sequences (
  sequence_id INTEGER PRIMARY KEY,
  genome_id INTEGER NOT NULL,        -- → reference_genomes(genome_id)
  sequence_name TEXT NOT NULL,       -- Chromosome/contig name
  sequence_length INTEGER NOT NULL,
  sequence BLOB,                     -- Actual sequence data
  FOREIGN KEY (genome_id) REFERENCES reference_genomes (genome_id)
)
```

**When populated:** `import_reference_genome()` splits FASTA into individual sequences

**How used:** `extract_flanking_sequences()` retrieves sequence data for coordinates

**Status:** ⚠️ **Legacy** - May be simplified in streamlined workflow

---

### 6. flanking_sequences ✅
**Purpose:** Extracted DNA sequences around variant positions

**Schema:**
```sql
CREATE TABLE flanking_sequences (
  flanking_id INTEGER PRIMARY KEY,
  vcf_id INTEGER NOT NULL,           -- → vcf_data(vcf_id)
  sequence_id INTEGER NOT NULL,      -- → reference_sequences(sequence_id)
  flank_size INTEGER NOT NULL,       -- bp on each side
  start_position INTEGER NOT NULL,
  end_position INTEGER NOT NULL,
  sequence TEXT NOT NULL,            -- Extracted sequence
  seq_type TEXT NOT NULL DEFAULT 'raw', -- 'raw' or 'orf'
  seq_length INTEGER,
  FOREIGN KEY (vcf_id) REFERENCES vcf_data (vcf_id),
  FOREIGN KEY (sequence_id) REFERENCES reference_sequences (sequence_id)
)
```

**When populated:** `extract_flanking_sequences()` extracts sequences around variants

**How used:** `run_blast()` uses sequences as BLAST queries

**Status:** ✅ **Keep** - Core functionality for sequence extraction

---

### 7. blast_parameters 🚫
**Purpose:** Parameters for BLAST runs (database, e-value, etc.)

**Schema:**
```sql
CREATE TABLE blast_parameters (
  blast_param_id INTEGER PRIMARY KEY,
  blast_type TEXT NOT NULL,          -- "blastp", "blastx", etc.
  db_name TEXT NOT NULL,
  db_path TEXT NOT NULL,
  e_value REAL NOT NULL,
  max_hits INTEGER NOT NULL,
  execution_date TEXT NOT NULL
)
```

**When populated:** `run_blast()` stores BLAST parameters

**How used:** Referenced by `blast_results`, `blast_database_metadata`

**Status:** 🚫 **REMOVE** - "One database = one analysis" means one BLAST run, don't need parameter tracking

---

### 8. blast_results ✅
**Purpose:** BLAST search results linking sequences to proteins

**Schema:**
```sql
CREATE TABLE blast_results (
  blast_result_id INTEGER PRIMARY KEY,
  blast_param_id INTEGER NOT NULL,  -- 🚫 REMOVE THIS
  flanking_id INTEGER NOT NULL,     -- → flanking_sequences(flanking_id)
  hit_accession TEXT NOT NULL,      -- UniProt accession
  hit_description TEXT,
  percent_identity REAL,
  alignment_length INTEGER,
  mismatches INTEGER,
  gap_openings INTEGER,
  query_start INTEGER,
  query_end INTEGER,
  subject_start INTEGER,
  subject_end INTEGER,
  e_value REAL,
  bit_score REAL,
  FOREIGN KEY (flanking_id) REFERENCES flanking_sequences (flanking_id)
)
```

**When populated:** `run_blast()` stores BLAST search results

**How used:** `retrieve_annotations()` uses hit_accession to get protein annotations

**Status:** ✅ **Keep** - Core functionality, but remove `blast_param_id` column and foreign key

---

### 9. blast_database_metadata 🚫
**Purpose:** Metadata about BLAST databases used

**Schema:**
```sql
CREATE TABLE blast_database_metadata (
  metadata_id INTEGER PRIMARY KEY,
  blast_param_id INTEGER NOT NULL,  -- → blast_parameters(blast_param_id)
  db_path TEXT NOT NULL,
  db_name TEXT NOT NULL,
  db_full_path TEXT NOT NULL,
  db_title TEXT,
  num_sequences INTEGER,
  total_length INTEGER,
  db_date TEXT,
  db_version TEXT,
  longest_sequence INTEGER,
  extraction_date TEXT NOT NULL,
  raw_output TEXT,
  FOREIGN KEY (blast_param_id) REFERENCES blast_parameters (blast_param_id)
)
```

**When populated:** `run_blast()` extracts database metadata

**How used:** Informational only

**Status:** 🚫 **REMOVE** - Unnecessary metadata in simplified workflow

---

### 10. annotations ✅
**Purpose:** UniProt protein annotation data

**Schema:**
```sql
CREATE TABLE annotations (
  annotation_id INTEGER PRIMARY KEY,
  blast_result_id INTEGER NOT NULL, -- → blast_results(blast_result_id)
  uniprot_accession TEXT NOT NULL,
  entry_name TEXT,
  gene_names TEXT,
  retrieval_date TEXT NOT NULL,
  FOREIGN KEY (blast_result_id) REFERENCES blast_results (blast_result_id)
)
```

**When populated:** `retrieve_annotations()` gets UniProt data for BLAST hits

**How used:** 
- Parent table for GO, KEGG, Pfam, InterPro, eggNOG annotations
- `ora()` functions query via annotations to get functional terms

**Status:** ✅ **Keep** - Core functionality

---

### 11. go_terms ✅
**Purpose:** Gene Ontology term annotations

**Schema:**
```sql
CREATE TABLE go_terms (
  go_term_id INTEGER PRIMARY KEY,
  annotation_id INTEGER NOT NULL,   -- → annotations(annotation_id)
  go_id TEXT NOT NULL,             -- "GO:0008150"
  go_term TEXT NOT NULL,           -- "biological_process"
  go_category TEXT NOT NULL,       -- "P", "F", "C"
  go_evidence TEXT,                -- Evidence code
  FOREIGN KEY (annotation_id) REFERENCES annotations (annotation_id)
)
```

**When populated:** `retrieve_annotations()` extracts GO terms from UniProt

**How used:** `ora()` functions perform GO enrichment analysis

**Status:** ✅ **Keep** - Core functionality

---

### 12. kegg_references ✅
**Purpose:** KEGG pathway annotations

**Schema:**
```sql
CREATE TABLE kegg_references (
  kegg_ref_id INTEGER PRIMARY KEY,
  annotation_id INTEGER NOT NULL,   -- → annotations(annotation_id)
  kegg_id TEXT NOT NULL,           -- "ko00010"
  pathway_name TEXT,               -- "Glycolysis / Gluconeogenesis"
  FOREIGN KEY (annotation_id) REFERENCES annotations (annotation_id)
)
```

**When populated:** `retrieve_annotations()` extracts KEGG pathways from UniProt

**How used:** `ora()` functions perform KEGG pathway enrichment

**Status:** ✅ **Keep** - Core functionality

---

### 13. pfam_domains ✅
**Purpose:** Pfam protein domain annotations

**Schema:**
```sql
CREATE TABLE pfam_domains (
  pfam_domain_id INTEGER PRIMARY KEY,
  annotation_id INTEGER NOT NULL,  -- → annotations(annotation_id)
  pfam_id TEXT NOT NULL,          -- "PF00001"
  domain_name TEXT,               -- "7 transmembrane receptor"
  match_status TEXT,              -- Quality of match
  FOREIGN KEY (annotation_id) REFERENCES annotations (annotation_id)
)
```

**When populated:** `retrieve_annotations()` extracts Pfam domains from UniProt

**How used:** Future Pfam-based enrichment analysis

**Status:** ✅ **Keep** - Extended functionality

---

### 14. interpro_families ✅
**Purpose:** InterPro protein family annotations

**Schema:**
```sql
CREATE TABLE interpro_families (
  interpro_family_id INTEGER PRIMARY KEY,
  annotation_id INTEGER NOT NULL, -- → annotations(annotation_id)
  interpro_id TEXT NOT NULL,     -- "IPR000001"
  family_name TEXT,              -- "Kringle"
  FOREIGN KEY (annotation_id) REFERENCES annotations (annotation_id)
)
```

**When populated:** `retrieve_annotations()` extracts InterPro families from UniProt

**How used:** Future InterPro-based enrichment analysis

**Status:** ✅ **Keep** - Extended functionality

---

### 15. eggnog_categories ✅
**Purpose:** eggNOG functional category annotations

**Schema:**
```sql
CREATE TABLE eggnog_categories (
  eggnog_category_id INTEGER PRIMARY KEY,
  annotation_id INTEGER NOT NULL, -- → annotations(annotation_id)
  eggnog_id TEXT NOT NULL,       -- "COG0001"
  taxonomic_scope TEXT,          -- Taxonomic level
  FOREIGN KEY (annotation_id) REFERENCES annotations (annotation_id)
)
```

**When populated:** `retrieve_annotations()` extracts eggNOG categories from UniProt

**How used:** Future eggNOG-based enrichment analysis

**Status:** ✅ **Keep** - Extended functionality

---

### 16. analysis_reports ⚠️
**Purpose:** Metadata about generated analysis reports

**Schema:**
```sql
CREATE TABLE analysis_reports (
  report_id INTEGER PRIMARY KEY,
  report_path TEXT NOT NULL,
  format TEXT NOT NULL,           -- "html", "pdf", etc.
  template TEXT NOT NULL,
  created_date TEXT NOT NULL,
  last_updated TEXT NOT NULL
)
```

**When populated:** Report generation functions

**How used:** Tracking generated reports

**Status:** ⚠️ **Review** - May not be needed in simplified workflow

---

### 17. ora_analyses ✅*
**Purpose:** Metadata for Over-Representation Analysis runs

**Schema:**
```sql
CREATE TABLE ora_analyses (
  analysis_id INTEGER PRIMARY KEY,
  blast_param_id INTEGER,          -- 🚫 REMOVE THIS
  annotation_type TEXT NOT NULL,   -- "GO", "KEGG"
  term_type TEXT NOT NULL,         -- "BP", "MF", "CC", "PATHWAY"
  analysis_date TEXT NOT NULL,
  total_foreground_genes INTEGER,
  total_background_genes INTEGER,
  analysis_parameters TEXT,        -- JSON parameters
  enrichment_method TEXT DEFAULT 'clusterprofiler',
  FOREIGN KEY (blast_param_id) REFERENCES blast_parameters (blast_param_id) -- 🚫 REMOVE
)
```

**When populated:** `store_ora_results()` during `ora()` analysis

**How used:** Referenced by `ora_results` table

**Status:** ✅ **Keep** - Core functionality, but remove `blast_param_id` column and foreign key

---

### 18. ora_results ✅
**Purpose:** Individual enrichment analysis results

**Schema:**
```sql
CREATE TABLE ora_results (
  result_id INTEGER PRIMARY KEY,
  analysis_id INTEGER NOT NULL,    -- → ora_analyses(analysis_id)
  term_id TEXT NOT NULL,          -- "GO:0008150", "ko00010"
  term_name TEXT NOT NULL,
  annotation_type TEXT NOT NULL,  -- "GO", "KEGG"
  term_type TEXT NOT NULL,        -- "BP", "MF", "CC", "PATHWAY"
  foreground_count INTEGER NOT NULL,
  background_count INTEGER NOT NULL,
  total_foreground INTEGER NOT NULL,
  total_background INTEGER NOT NULL,
  expected_count REAL,
  fold_enrichment REAL,
  p_value REAL NOT NULL,
  p_adjusted REAL,               -- FDR corrected
  significance_level TEXT,       -- "significant", "not_significant"
  gene_ratio TEXT,              -- clusterProfiler format
  bg_ratio TEXT,                -- clusterProfiler format
  qvalue REAL,                  -- clusterProfiler q-value
  gene_ids TEXT,                -- Comma-separated gene list
  FOREIGN KEY (analysis_id) REFERENCES ora_analyses (analysis_id)
)
```

**When populated:** `store_ora_results()` during `ora()` analysis

**How used:** Query functions retrieve enrichment results for visualization and analysis

**Status:** ✅ **Keep** - Core functionality

---

### 19. uniprot_cache ✅
**Purpose:** Cache UniProt API responses to avoid repeated requests

**Schema:**
```sql
CREATE TABLE uniprot_cache (
  cache_id INTEGER PRIMARY KEY,
  accession TEXT NOT NULL UNIQUE, -- UniProt accession
  response_json TEXT NOT NULL,    -- Full API response
  retrieval_date TEXT NOT NULL,
  UNIQUE(accession)
)
```

**When populated:** `retrieve_annotations()` caches UniProt API responses

**How used:** `retrieve_annotations()` checks cache before making API calls

**Status:** ✅ **Keep** - Performance optimization

---

### 20. locus_statistics ✅
**Purpose:** Statistical values (p-values, q-values) for genomic loci

**Schema:**
```sql
CREATE TABLE locus_statistics (
  statistic_id INTEGER PRIMARY KEY,
  chromosome TEXT NOT NULL,
  position INTEGER NOT NULL,
  statistic REAL NOT NULL,        -- p-value, q-value, etc.
  created_date TEXT NOT NULL,
  UNIQUE(chromosome, position)    -- Prevents duplicates
)
```

**When populated:** `define_locus_statistics()` in streamlined workflow

**How used:** 
- `define_locus_statistics()` with threshold creates candidates
- Future: extract p-values for ReviGO export

**Status:** ✅ **Keep** - Core streamlined workflow functionality

---

### 21. candidate_loci ✅
**Purpose:** Defined candidate loci for enrichment analysis

**Schema:**
```sql
CREATE TABLE candidate_loci (
  candidate_id INTEGER PRIMARY KEY,
  chromosome TEXT NOT NULL,
  position INTEGER NOT NULL,
  method TEXT NOT NULL,           -- "threshold", "manual", etc.
  threshold REAL,                 -- Threshold used (if applicable)
  created_date TEXT NOT NULL,
  UNIQUE(chromosome, position)    -- Prevents duplicates
)
```

**When populated:** 
- `define_locus_statistics()` with threshold
- `define_candidate_loci()` manual definition

**How used:** `ora()` functions use stored candidates for enrichment analysis

**Status:** ✅ **Keep** - Core streamlined workflow functionality

---

## Legacy Cleanup Summary

### 🚫 **Tables to Remove Completely:**
1. `blast_parameters` - One database = one analysis, no parameter tracking needed
2. `blast_database_metadata` - Unnecessary metadata

### 🚫 **Columns to Remove:**
1. `blast_results.blast_param_id` + foreign key constraint
2. `ora_analyses.blast_param_id` + foreign key constraint

### ⚠️ **Tables to Review:**
1. `input_files` - May not be needed in streamlined workflow
2. `vcf_data` - May not be needed if users provide coordinates directly
3. `reference_genomes` + `reference_sequences` - May be simplified
4. `analysis_reports` - May not be needed

### ✅ **Core Tables (Keep):**
- `metadata` - Database validation
- `flanking_sequences` - Sequence extraction
- `blast_results` - BLAST results (minus blast_param_id)
- `annotations` + all annotation child tables - Functional annotations
- `ora_analyses` + `ora_results` - Enrichment analysis (minus blast_param_id)
- `uniprot_cache` - Performance
- `locus_statistics` + `candidate_loci` - Streamlined workflow