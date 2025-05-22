#' Database schema definition for funseqR
#'
#' This file contains functions to create and manage the schema of the funseqR database.
#'

#' Create the funseqR database schema
#'
#' This function creates the necessary tables and indexes for the funseqR database.
#'
#' @param con A database connection object.
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE.
#'
#' @return Invisible NULL.
#'
#' @importFrom DBI dbExecute
#' @export
create_funseq_schema <- function(con, verbose = TRUE) {
  # Create tables
  if (verbose) message("Creating metadata table...")
  DBI::dbExecute(con, "
    CREATE TABLE metadata (
      key TEXT PRIMARY KEY,
      value TEXT NOT NULL
    )
  ")

  if (verbose) message("Creating projects table...")
  DBI::dbExecute(con, "
    CREATE TABLE projects (
      project_id INTEGER PRIMARY KEY,
      project_name TEXT NOT NULL,
      description TEXT,
      creation_date TEXT NOT NULL,
      last_modified TEXT NOT NULL
    )
  ")

  if (verbose) message("Creating input_files table...")
  DBI::dbExecute(con, "
    CREATE TABLE input_files (
      file_id INTEGER PRIMARY KEY,
      project_id INTEGER NOT NULL,
      file_type TEXT NOT NULL,
      file_name TEXT NOT NULL,
      file_path TEXT NOT NULL,
      file_hash TEXT,
      import_date TEXT NOT NULL,
      FOREIGN KEY (project_id) REFERENCES projects (project_id)
    )
  ")

  if (verbose) message("Creating vcf_data table...")
  DBI::dbExecute(con, "
    CREATE TABLE vcf_data (
      vcf_id INTEGER PRIMARY KEY,
      file_id INTEGER NOT NULL,
      chromosome TEXT NOT NULL,
      position INTEGER NOT NULL,
      id TEXT,
      ref TEXT NOT NULL,
      alt TEXT NOT NULL,
      qual REAL,
      filter TEXT,
      info TEXT,
      format TEXT,
      sample_data TEXT,
      FOREIGN KEY (file_id) REFERENCES input_files (file_id)
    )
  ")

  if (verbose) message("Creating reference_genomes table...")
  DBI::dbExecute(con, "
    CREATE TABLE reference_genomes (
      genome_id INTEGER PRIMARY KEY,
      file_id INTEGER NOT NULL,
      genome_name TEXT NOT NULL,
      genome_build TEXT,
      FOREIGN KEY (file_id) REFERENCES input_files (file_id)
    )
  ")

  if (verbose) message("Creating reference_sequences table...")
  DBI::dbExecute(con, "
    CREATE TABLE reference_sequences (
      sequence_id INTEGER PRIMARY KEY,
      genome_id INTEGER NOT NULL,
      sequence_name TEXT NOT NULL,
      sequence_length INTEGER NOT NULL,
      sequence BLOB,
      FOREIGN KEY (genome_id) REFERENCES reference_genomes (genome_id)
    )
  ")

  if (verbose) message("Creating flanking_sequences table...")
  DBI::dbExecute(con, "
    CREATE TABLE flanking_sequences (
      flanking_id INTEGER PRIMARY KEY,
      vcf_id INTEGER NOT NULL,
      sequence_id INTEGER NOT NULL,
      flank_size INTEGER NOT NULL,
      start_position INTEGER NOT NULL,
      end_position INTEGER NOT NULL,
      sequence TEXT NOT NULL,
      FOREIGN KEY (vcf_id) REFERENCES vcf_data (vcf_id),
      FOREIGN KEY (sequence_id) REFERENCES reference_sequences (sequence_id)
    )
  ")

  if (verbose) message("Creating blast_parameters table...")
  DBI::dbExecute(con, "
    CREATE TABLE blast_parameters (
      blast_param_id INTEGER PRIMARY KEY,
      project_id INTEGER NOT NULL,
      blast_type TEXT NOT NULL,
      db_name TEXT NOT NULL,
      db_path TEXT NOT NULL,
      e_value REAL NOT NULL,
      max_hits INTEGER NOT NULL,
      execution_date TEXT NOT NULL,
      FOREIGN KEY (project_id) REFERENCES projects (project_id)
    )
  ")

  if (verbose) message("Creating blast_results table...")
  DBI::dbExecute(con, "
    CREATE TABLE blast_results (
      blast_result_id INTEGER PRIMARY KEY,
      blast_param_id INTEGER NOT NULL,
      flanking_id INTEGER NOT NULL,
      hit_accession TEXT NOT NULL,
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
      FOREIGN KEY (blast_param_id) REFERENCES blast_parameters (blast_param_id),
      FOREIGN KEY (flanking_id) REFERENCES flanking_sequences (flanking_id)
    )
  ")

  if (verbose) message("Creating blast_database_metadata table...")
  DBI::dbExecute(con, "
    CREATE TABLE blast_database_metadata (
      metadata_id INTEGER PRIMARY KEY,
      blast_param_id INTEGER NOT NULL,
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
  ")

  if (verbose) message("Creating annotations table...")
  DBI::dbExecute(con, "
    CREATE TABLE annotations (
      annotation_id INTEGER PRIMARY KEY,
      blast_result_id INTEGER NOT NULL,
      uniprot_accession TEXT NOT NULL,
      entry_name TEXT,
      gene_names TEXT,
      retrieval_date TEXT NOT NULL,
      FOREIGN KEY (blast_result_id) REFERENCES blast_results (blast_result_id)
    )
  ")

  if (verbose) message("Creating go_terms table...")
  DBI::dbExecute(con, "
    CREATE TABLE go_terms (
      go_term_id INTEGER PRIMARY KEY,
      annotation_id INTEGER NOT NULL,
      go_id TEXT NOT NULL,
      go_term TEXT NOT NULL,
      go_category TEXT NOT NULL,
      go_evidence TEXT,
      FOREIGN KEY (annotation_id) REFERENCES annotations (annotation_id)
    )
  ")

  if (verbose) message("Creating kegg_references table...")
  DBI::dbExecute(con, "
    CREATE TABLE kegg_references (
      kegg_ref_id INTEGER PRIMARY KEY,
      annotation_id INTEGER NOT NULL,
      kegg_id TEXT NOT NULL,
      pathway_name TEXT,
      FOREIGN KEY (annotation_id) REFERENCES annotations (annotation_id)
    )
  ")

  # Create indexes
  if (verbose) message("Creating indexes...")

  # VCF data indexes
  DBI::dbExecute(con, "CREATE INDEX idx_vcf_chrom_pos ON vcf_data (chromosome, position)")
  DBI::dbExecute(con, "CREATE INDEX idx_vcf_file_id ON vcf_data (file_id)")

  # Sequences indexes
  DBI::dbExecute(con, "CREATE INDEX idx_ref_seq_name ON reference_sequences (sequence_name)")
  DBI::dbExecute(con, "CREATE INDEX idx_ref_seq_genome ON reference_sequences (genome_id)")

  # Flanking sequences indexes
  DBI::dbExecute(con, "CREATE INDEX idx_flanking_vcf_id ON flanking_sequences (vcf_id)")
  DBI::dbExecute(con, "CREATE INDEX idx_flanking_seq_id ON flanking_sequences (sequence_id)")

  # BLAST results indexes
  DBI::dbExecute(con, "CREATE INDEX idx_blast_res_param ON blast_results (blast_param_id)")
  DBI::dbExecute(con, "CREATE INDEX idx_blast_res_flanking ON blast_results (flanking_id)")
  DBI::dbExecute(con, "CREATE INDEX idx_blast_res_accession ON blast_results (hit_accession)")

  # Annotation indexes
  DBI::dbExecute(con, "CREATE INDEX idx_anno_blast_result ON annotations (blast_result_id)")
  DBI::dbExecute(con, "CREATE INDEX idx_anno_uniprot ON annotations (uniprot_accession)")

  # GO terms and KEGG indexes
  DBI::dbExecute(con, "CREATE INDEX idx_go_annotation ON go_terms (annotation_id)")
  DBI::dbExecute(con, "CREATE INDEX idx_go_id ON go_terms (go_id)")
  DBI::dbExecute(con, "CREATE INDEX idx_go_category ON go_terms (go_category)")
  DBI::dbExecute(con, "CREATE INDEX idx_kegg_annotation ON kegg_references (annotation_id)")
  DBI::dbExecute(con, "CREATE INDEX idx_kegg_id ON kegg_references (kegg_id)")

  # BLAST database metadata indexes
  DBI::dbExecute(con, "CREATE INDEX idx_blast_db_meta_param ON blast_database_metadata (blast_param_id)")
  DBI::dbExecute(con, "CREATE INDEX idx_blast_db_meta_name ON blast_database_metadata (db_name)")
  DBI::dbExecute(con, "CREATE INDEX idx_blast_db_meta_date ON blast_database_metadata (extraction_date)")

  if (verbose) message("Schema creation complete.")

  return(invisible(NULL))
}

#' Upgrade the funseqR database schema
#'
#' This function upgrades the schema of an existing funseqR database to the current version.
#'
#' @param con A database connection object.
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE.
#'
#' @return Invisible NULL.
#'
#' @importFrom DBI dbGetQuery dbExecute
#' @export
upgrade_funseq_schema <- function(con, verbose = TRUE) {
  # Get current schema version
  if ("metadata" %in% DBI::dbListTables(con)) {
    metadata <- DBI::dbGetQuery(con, "SELECT key, value FROM metadata WHERE key = 'schema_version'")

    if (nrow(metadata) > 0) {
      current_version <- as.numeric(metadata$value)
    } else {
      current_version <- 1.0  # Assume version 1.0 if not specified
    }
  } else {
    stop("This does not appear to be a valid funseqR database.")
  }

  # Define the current schema version
  latest_version <- 1.0

  if (current_version >= latest_version) {
    if (verbose) message("Database schema is already at the latest version.")
    return(invisible(NULL))
  }

  # Upgrade path for each version
  if (current_version < 1.0) {
    # Upgrade to version 1.0
    if (verbose) message("Upgrading schema to version 1.0...")

    # Add version tracking to metadata
    if (DBI::dbGetQuery(con, "SELECT COUNT(*) FROM metadata WHERE key = 'schema_version'")[[1]] == 0) {
      DBI::dbExecute(con, "INSERT INTO metadata (key, value) VALUES ('schema_version', '1.0')")
    } else {
      DBI::dbExecute(con, "UPDATE metadata SET value = '1.0' WHERE key = 'schema_version'")
    }
  }

  if (verbose) message("Schema upgrade complete.")

  return(invisible(NULL))
}

#' Check the schema version of a funseqR database
#'
#' @param con A database connection object.
#'
#' @return A numeric value representing the schema version.
#'
#' @importFrom DBI dbGetQuery
#' @export
get_schema_version <- function(con) {
  if ("metadata" %in% DBI::dbListTables(con)) {
    metadata <- DBI::dbGetQuery(con, "SELECT value FROM metadata WHERE key = 'schema_version'")

    if (nrow(metadata) > 0) {
      return(as.numeric(metadata$value))
    } else {
      return(1.0)  # Assume version 1.0 if not specified
    }
  } else {
    stop("This does not appear to be a valid funseqR database.")
  }
}
