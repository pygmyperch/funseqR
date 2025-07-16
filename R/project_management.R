#' Project Management Functions for funseqR
#'
#' Functions for managing funseqR projects, including data export and project duplication.
#'

#' Export data from existing project to create new project
#'
#' Creates a new funseqR database by selectively copying data from an existing project.
#' This enables creating new analysis projects without repeating time-consuming 
#' computational steps (VCF import, reference import, flanking sequences, BLAST, annotations).
#'
#' @param source_con Database connection object for the source project
#' @param target_db_path Character string specifying the path for the new database file
#' @param export_vcf Logical. Export VCF data and file registrations. Default is TRUE
#' @param export_reference Logical. Export reference genome data. Default is TRUE
#' @param export_chromosomes Logical. Export chromosome definitions from define_chromosomes(). Default is TRUE
#' @param export_candidates Logical. Export candidate loci and associated locus statistics (if available). Default is FALSE
#' @param export_flanking Logical. Export flanking sequences. Default is FALSE
#' @param export_blast Logical. Export BLAST parameters, results, and database metadata. Default is FALSE
#' @param export_annotations Logical. Export all annotation data (GO, KEGG, Pfam, InterPro, eggNOG). Default is FALSE
#' @param force Logical. If TRUE, overwrite existing database file. Default is FALSE
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return Database connection object to the newly created project database
#'
#' @details
#' This function implements the "one database = one analysis project" philosophy by enabling
#' selective data export between projects. It's designed to avoid repeating time-consuming
#' computational steps while maintaining data integrity and foreign key relationships.
#'
#' \strong{Export Categories:}
#'
#' \strong{Essential Data (fast to export):}
#' \itemize{
#'   \item \code{export_vcf}: VCF data, coordinates, sample information
#'   \item \code{export_reference}: Reference genome sequences and metadata
#'   \item \code{export_chromosomes}: Main chromosome definitions for visualization
#' }
#'
#' \strong{Analysis Data:}
#' \itemize{
#'   \item \code{export_candidates}: Candidate loci coordinates + underlying locus statistics (when available)
#' }
#'
#' \strong{Computational Data (time-intensive to recreate):}
#' \itemize{
#'   \item \code{export_flanking}: Flanking sequences around variants (slow: depends on sequence extraction)
#'   \item \code{export_blast}: BLAST search results (very slow: depends on BLAST runtime)
#'   \item \code{export_annotations}: Functional annotations (slow: depends on annotation processing)
#' }
#'
#' \strong{Smart Dependencies:}
#' When \code{export_candidates = TRUE}, the function automatically includes locus statistics
#' if they exist, ensuring the new project can fully reproduce candidate analysis workflows.
#'
#' \strong{Use Cases:}
#' \itemize{
#'   \item \strong{Analysis branching}: Create variants of a project for different analysis approaches
#'   \item \strong{Collaboration}: Share base data without computational results
#'   \item \strong{Method comparison}: Test different annotation or enrichment approaches
#'   \item \strong{Data archiving}: Preserve essential data while experimenting
#' }
#'
#' @examples
#' \dontrun{
#' # Connect to existing project
#' source_con <- connect_funseq_db("original_project.db")
#' 
#' # Minimal export - just essential data (fast)
#' new_con <- export_project_data(source_con, "new_analysis.db")
#' 
#' # Include candidate analysis capability
#' new_con <- export_project_data(
#'   source_con, 
#'   "candidate_analysis.db", 
#'   export_candidates = TRUE  # Auto-includes locus_statistics if available
#' )
#' 
#' # Complete computational export (all data)
#' new_con <- export_project_data(
#'   source_con,
#'   "complete_copy.db",
#'   export_candidates = TRUE,
#'   export_flanking = TRUE,
#'   export_blast = TRUE,
#'   export_annotations = TRUE
#' )
#' 
#' # Custom export - keep sequences, re-run BLAST and annotations
#' new_con <- export_project_data(
#'   source_con,
#'   "reanalysis.db", 
#'   export_candidates = TRUE,
#'   export_flanking = TRUE,   # Keep sequences
#'   export_blast = FALSE,     # Re-run BLAST with different parameters
#'   export_annotations = FALSE # Re-annotate with updated databases
#' )
#' }
#'
#' @export
export_project_data <- function(source_con, target_db_path,
                                export_vcf = TRUE,
                                export_reference = TRUE,
                                export_chromosomes = TRUE,
                                export_candidates = FALSE,
                                export_flanking = FALSE,
                                export_blast = FALSE,
                                export_annotations = FALSE,
                                force = FALSE,
                                verbose = TRUE) {
  
  if (verbose) message("=== Exporting Project Data ===")
  if (verbose) message("Source: ", DBI::dbGetInfo(source_con)$dbname)
  if (verbose) message("Target: ", target_db_path)
  
  # Validate source connection
  if (!DBI::dbIsValid(source_con)) {
    stop("Invalid source database connection")
  }
  
  # Create target database
  if (verbose) message("\n--- Creating target database ---")
  target_con <- create_funseq_db(target_db_path, force = force, verbose = verbose)
  
  tryCatch({
    # Track ID mappings for foreign key updates
    id_mappings <- list()
    
    # Export essential data
    if (export_vcf) {
      if (verbose) message("\n--- Exporting VCF data ---")
      id_mappings <- .export_vcf_data(source_con, target_con, id_mappings, verbose)
    }
    
    if (export_reference) {
      if (verbose) message("\n--- Exporting reference genome data ---")
      id_mappings <- .export_reference_data(source_con, target_con, id_mappings, verbose)
    }
    
    if (export_chromosomes) {
      if (verbose) message("\n--- Exporting chromosome definitions ---")
      .export_chromosome_metadata(source_con, target_con, verbose)
    }
    
    # Export candidate data with smart dependencies
    if (export_candidates) {
      if (verbose) message("\n--- Exporting candidate analysis data ---")
      .export_candidate_data(source_con, target_con, verbose)
    }
    
    # Export computational data
    if (export_flanking) {
      if (verbose) message("\n--- Exporting flanking sequences ---")
      .export_flanking_data(source_con, target_con, id_mappings, verbose)
    }
    
    if (export_blast) {
      if (verbose) message("\n--- Exporting BLAST data ---")
      id_mappings <- .export_blast_data(source_con, target_con, id_mappings, verbose)
    }
    
    if (export_annotations) {
      if (verbose) message("\n--- Exporting annotation data ---")
      .export_annotation_data(source_con, target_con, id_mappings, verbose)
    }
    
    if (verbose) {
      message("\n=== Export completed successfully ===")
      message("New project database: ", target_db_path)
      
      # Show summary of exported data
      .show_export_summary(target_con, verbose)
    }
    
    return(target_con)
    
  }, error = function(e) {
    # Clean up on error
    if (DBI::dbIsValid(target_con)) {
      DBI::dbDisconnect(target_con)
    }
    if (file.exists(target_db_path)) {
      file.remove(target_db_path)
      if (verbose) message("Removed incomplete database file due to error")
    }
    stop("Export failed: ", e$message)
  })
}

#' Export VCF data and file registrations
#' @keywords internal
.export_vcf_data <- function(source_con, target_con, id_mappings, verbose) {
  
  # Get VCF files from source
  vcf_files <- DBI::dbGetQuery(source_con, 
    "SELECT * FROM input_files WHERE file_type = 'vcf' ORDER BY file_id")
  
  if (nrow(vcf_files) == 0) {
    if (verbose) message("  - No VCF files found in source database")
    return(id_mappings)
  }
  
  if (verbose) message("  - Exporting ", nrow(vcf_files), " VCF file(s)")
  
  # Initialize file ID mapping
  id_mappings$file_ids <- data.frame(
    old_id = integer(0),
    new_id = integer(0)
  )
  
  # Export VCF files and data
  for (i in 1:nrow(vcf_files)) {
    file_record <- vcf_files[i, ]
    old_file_id <- file_record$file_id
    
    # Insert file record (auto-generates new file_id)
    DBI::dbExecute(target_con, "
      INSERT INTO input_files (file_type, file_name, file_path, file_hash, import_date)
      VALUES (?, ?, ?, ?, ?)
    ", list(file_record$file_type, file_record$file_name, file_record$file_path,
            file_record$file_hash, file_record$import_date))
    
    # Get new file_id
    new_file_id <- DBI::dbGetQuery(target_con, "SELECT last_insert_rowid() as id")$id
    
    # Store mapping
    id_mappings$file_ids <- rbind(id_mappings$file_ids, 
                                  data.frame(old_id = old_file_id, new_id = new_file_id))
    
    # Export VCF data for this file
    vcf_data <- DBI::dbGetQuery(source_con,
      "SELECT * FROM vcf_data WHERE file_id = ? ORDER BY vcf_id", 
      list(old_file_id))
    
    if (nrow(vcf_data) > 0) {
      if (verbose) message("    - Exporting ", nrow(vcf_data), " VCF entries for file: ", file_record$file_name)
      
      # Update file_id and insert VCF data
      vcf_data$file_id <- new_file_id
      
      # Batch insert VCF data
      .batch_insert_vcf_data(target_con, vcf_data)
    }
  }
  
  return(id_mappings)
}

#' Batch insert VCF data
#' @keywords internal
.batch_insert_vcf_data <- function(con, vcf_data) {
  
  # Insert in batches for efficiency
  batch_size <- 1000
  num_batches <- ceiling(nrow(vcf_data) / batch_size)
  
  for (batch in 1:num_batches) {
    start_idx <- (batch - 1) * batch_size + 1
    end_idx <- min(batch * batch_size, nrow(vcf_data))
    batch_data <- vcf_data[start_idx:end_idx, ]
    
    # Create parameterized query for batch
    placeholders <- paste0("(", paste(rep("?", 11), collapse = ", "), ")")
    values_clause <- paste(rep(placeholders, nrow(batch_data)), collapse = ", ")
    query <- paste0(
      "INSERT INTO vcf_data (file_id, chromosome, position, id, ref, alt, qual, filter, info, format, sample_data) VALUES ",
      values_clause
    )
    
    # Flatten parameters for batch insert
    params <- list()
    for (i in 1:nrow(batch_data)) {
      params <- c(params, list(
        batch_data$file_id[i],
        batch_data$chromosome[i],
        batch_data$position[i],
        batch_data$id[i],
        batch_data$ref[i],
        batch_data$alt[i],
        batch_data$qual[i],
        batch_data$filter[i],
        batch_data$info[i],
        batch_data$format[i],
        batch_data$sample_data[i]
      ))
    }
    
    DBI::dbExecute(con, query, params)
  }
}

#' Export reference genome data
#' @keywords internal
.export_reference_data <- function(source_con, target_con, id_mappings, verbose) {
  
  # Get reference files from source
  ref_files <- DBI::dbGetQuery(source_con, 
    "SELECT * FROM input_files WHERE file_type = 'fasta' ORDER BY file_id")
  
  if (nrow(ref_files) == 0) {
    if (verbose) message("  - No reference files found in source database")
    return(id_mappings)
  }
  
  if (verbose) message("  - Exporting ", nrow(ref_files), " reference file(s)")
  
  # Export reference files and genomes
  for (i in 1:nrow(ref_files)) {
    file_record <- ref_files[i, ]
    old_file_id <- file_record$file_id
    
    # Insert file record
    DBI::dbExecute(target_con, "
      INSERT INTO input_files (file_type, file_name, file_path, file_hash, import_date)
      VALUES (?, ?, ?, ?, ?)
    ", list(file_record$file_type, file_record$file_name, file_record$file_path,
            file_record$file_hash, file_record$import_date))
    
    # Get new file_id
    new_file_id <- DBI::dbGetQuery(target_con, "SELECT last_insert_rowid() as id")$id
    
    # Update file ID mapping
    id_mappings$file_ids <- rbind(id_mappings$file_ids, 
                                  data.frame(old_id = old_file_id, new_id = new_file_id))
    
    # Export reference genomes for this file
    ref_genomes <- DBI::dbGetQuery(source_con,
      "SELECT * FROM reference_genomes WHERE file_id = ?", 
      list(old_file_id))
    
    if (nrow(ref_genomes) > 0) {
      # Initialize genome ID mapping if not exists
      if (!"genome_ids" %in% names(id_mappings)) {
        id_mappings$genome_ids <- data.frame(old_id = integer(0), new_id = integer(0))
      }
      
      for (j in 1:nrow(ref_genomes)) {
        genome_record <- ref_genomes[j, ]
        old_genome_id <- genome_record$genome_id
        
        # Insert genome record
        DBI::dbExecute(target_con, "
          INSERT INTO reference_genomes (file_id, genome_name, genome_build)
          VALUES (?, ?, ?)
        ", list(new_file_id, genome_record$genome_name, genome_record$genome_build))
        
        # Get new genome_id
        new_genome_id <- DBI::dbGetQuery(target_con, "SELECT last_insert_rowid() as id")$id
        
        # Store genome mapping
        id_mappings$genome_ids <- rbind(id_mappings$genome_ids,
                                        data.frame(old_id = old_genome_id, new_id = new_genome_id))
        
        # Export reference sequences for this genome
        ref_sequences <- DBI::dbGetQuery(source_con,
          "SELECT * FROM reference_sequences WHERE genome_id = ?", 
          list(old_genome_id))
        
        if (nrow(ref_sequences) > 0) {
          if (verbose) message("    - Exporting ", nrow(ref_sequences), " reference sequences for genome: ", genome_record$genome_name)
          
          # Update genome_id and insert sequences
          ref_sequences$genome_id <- new_genome_id
          
          # Batch insert sequences
          .batch_insert_reference_sequences(target_con, ref_sequences)
        }
      }
    }
  }
  
  return(id_mappings)
}

#' Batch insert reference sequences
#' @keywords internal
.batch_insert_reference_sequences <- function(con, ref_sequences) {
  
  # Insert in batches for efficiency
  batch_size <- 100
  num_batches <- ceiling(nrow(ref_sequences) / batch_size)
  
  for (batch in 1:num_batches) {
    start_idx <- (batch - 1) * batch_size + 1
    end_idx <- min(batch * batch_size, nrow(ref_sequences))
    batch_data <- ref_sequences[start_idx:end_idx, ]
    
    # Create parameterized query for batch
    placeholders <- paste0("(", paste(rep("?", 4), collapse = ", "), ")")
    values_clause <- paste(rep(placeholders, nrow(batch_data)), collapse = ", ")
    query <- paste0(
      "INSERT INTO reference_sequences (genome_id, sequence_name, sequence_length, sequence) VALUES ",
      values_clause
    )
    
    # Flatten parameters for batch insert
    params <- list()
    for (i in 1:nrow(batch_data)) {
      params <- c(params, list(
        batch_data$genome_id[i],
        batch_data$sequence_name[i],
        batch_data$sequence_length[i],
        batch_data$sequence[i]
      ))
    }
    
    DBI::dbExecute(con, query, params)
  }
}

#' Export chromosome metadata
#' @keywords internal
.export_chromosome_metadata <- function(source_con, target_con, verbose) {
  
  # Get chromosome definitions from metadata
  chr_metadata <- DBI::dbGetQuery(source_con,
    "SELECT * FROM metadata WHERE key = 'main_chromosomes'")
  
  if (nrow(chr_metadata) > 0) {
    if (verbose) message("  - Exporting chromosome definitions")
    
    # Insert chromosome metadata
    DBI::dbExecute(target_con, "
      INSERT INTO metadata (key, value) VALUES (?, ?)
    ", list(chr_metadata$key[1], chr_metadata$value[1]))
    
  } else {
    if (verbose) message("  - No chromosome definitions found in source database")
  }
}

#' Export candidate data with smart dependencies
#' @keywords internal
.export_candidate_data <- function(source_con, target_con, verbose) {
  
  # Check for candidate loci
  candidates <- DBI::dbGetQuery(source_con, "SELECT * FROM candidate_loci ORDER BY candidate_id")
  
  if (nrow(candidates) > 0) {
    if (verbose) message("  - Exporting ", nrow(candidates), " candidate loci")
    
    # Export candidate loci
    .batch_insert_candidate_loci(target_con, candidates)
    
    # Check for locus statistics (smart dependency)
    statistics <- DBI::dbGetQuery(source_con, "SELECT COUNT(*) as count FROM locus_statistics")$count
    
    if (statistics > 0) {
      if (verbose) message("  - Exporting ", statistics, " locus statistics (candidate dependency)")
      
      # Export all locus statistics
      locus_stats <- DBI::dbGetQuery(source_con, "SELECT * FROM locus_statistics ORDER BY statistic_id")
      .batch_insert_locus_statistics(target_con, locus_stats)
    } else {
      if (verbose) message("  - No locus statistics found (candidates may have been defined via coordinates)")
    }
    
  } else {
    if (verbose) message("  - No candidate loci found in source database")
  }
}

#' Batch insert candidate loci
#' @keywords internal
.batch_insert_candidate_loci <- function(con, candidates) {
  
  # Insert in batches
  batch_size <- 1000
  num_batches <- ceiling(nrow(candidates) / batch_size)
  
  for (batch in 1:num_batches) {
    start_idx <- (batch - 1) * batch_size + 1
    end_idx <- min(batch * batch_size, nrow(candidates))
    batch_data <- candidates[start_idx:end_idx, ]
    
    # Create parameterized query for batch
    placeholders <- paste0("(", paste(rep("?", 5), collapse = ", "), ")")
    values_clause <- paste(rep(placeholders, nrow(batch_data)), collapse = ", ")
    query <- paste0(
      "INSERT INTO candidate_loci (chromosome, position, method, threshold, created_date) VALUES ",
      values_clause
    )
    
    # Flatten parameters for batch insert
    params <- list()
    for (i in 1:nrow(batch_data)) {
      params <- c(params, list(
        batch_data$chromosome[i],
        batch_data$position[i],
        batch_data$method[i],
        batch_data$threshold[i],
        batch_data$created_date[i]
      ))
    }
    
    DBI::dbExecute(con, query, params)
  }
}

#' Batch insert locus statistics
#' @keywords internal
.batch_insert_locus_statistics <- function(con, locus_stats) {
  
  # Insert in batches
  batch_size <- 1000
  num_batches <- ceiling(nrow(locus_stats) / batch_size)
  
  for (batch in 1:num_batches) {
    start_idx <- (batch - 1) * batch_size + 1
    end_idx <- min(batch * batch_size, nrow(locus_stats))
    batch_data <- locus_stats[start_idx:end_idx, ]
    
    # Create parameterized query for batch
    placeholders <- paste0("(", paste(rep("?", 4), collapse = ", "), ")")
    values_clause <- paste(rep(placeholders, nrow(batch_data)), collapse = ", ")
    query <- paste0(
      "INSERT INTO locus_statistics (chromosome, position, statistic, created_date) VALUES ",
      values_clause
    )
    
    # Flatten parameters for batch insert
    params <- list()
    for (i in 1:nrow(batch_data)) {
      params <- c(params, list(
        batch_data$chromosome[i],
        batch_data$position[i],
        batch_data$statistic[i],
        batch_data$created_date[i]
      ))
    }
    
    DBI::dbExecute(con, query, params)
  }
}

#' Export flanking sequences data
#' @keywords internal
.export_flanking_data <- function(source_con, target_con, id_mappings, verbose) {
  
  # Check for flanking sequences
  flanking_count <- DBI::dbGetQuery(source_con, "SELECT COUNT(*) as count FROM flanking_sequences")$count
  
  if (flanking_count == 0) {
    if (verbose) message("  - No flanking sequences found in source database")
    return(NULL)
  }
  
  if (verbose) message("  - Exporting ", flanking_count, " flanking sequences")
  
  # Get all flanking sequences with VCF mapping
  flanking_data <- DBI::dbGetQuery(source_con, "
    SELECT fs.*, vd.file_id as old_file_id
    FROM flanking_sequences fs
    JOIN vcf_data vd ON fs.vcf_id = vd.vcf_id
    ORDER BY fs.flanking_id
  ")
  
  if (nrow(flanking_data) > 0) {
    # Update vcf_id based on new VCF data
    flanking_export <- data.frame()
    
    for (i in 1:nrow(flanking_data)) {
      fs_record <- flanking_data[i, ]
      
      # Find corresponding VCF entry in target database
      old_file_id <- fs_record$old_file_id
      new_file_id <- id_mappings$file_ids$new_id[id_mappings$file_ids$old_id == old_file_id]
      
      if (length(new_file_id) > 0) {
        # Find new vcf_id
        new_vcf <- DBI::dbGetQuery(target_con, "
          SELECT vcf_id FROM vcf_data 
          WHERE file_id = ? AND chromosome = ? AND position = ?
        ", list(new_file_id[1], fs_record$chromosome, fs_record$position))
        
        if (nrow(new_vcf) > 0) {
          # Add to export data with updated vcf_id (use original column structure)
          export_record <- data.frame(
            vcf_id = new_vcf$vcf_id[1],
            chromosome = fs_record$chromosome,
            position = fs_record$position,
            start_position = fs_record$start_position,
            end_position = fs_record$end_position,
            raw_sequence = fs_record$raw_sequence,
            orf_nucleotide = fs_record$orf_nucleotide,
            orf_amino_acid = fs_record$orf_amino_acid,
            created_date = fs_record$created_date
          )
          
          flanking_export <- rbind(flanking_export, export_record)
        }
      }
    }
    
    if (nrow(flanking_export) > 0) {
      # Batch insert flanking sequences
      .batch_insert_flanking_sequences(target_con, flanking_export)
      if (verbose) message("    - Successfully exported ", nrow(flanking_export), " flanking sequences")
    }
  }
}

#' Batch insert flanking sequences
#' @keywords internal
.batch_insert_flanking_sequences <- function(con, flanking_data) {
  
  # Insert in batches
  batch_size <- 500
  num_batches <- ceiling(nrow(flanking_data) / batch_size)
  
  for (batch in 1:num_batches) {
    start_idx <- (batch - 1) * batch_size + 1
    end_idx <- min(batch * batch_size, nrow(flanking_data))
    batch_data <- flanking_data[start_idx:end_idx, ]
    
    # Create parameterized query for batch
    placeholders <- paste0("(", paste(rep("?", 9), collapse = ", "), ")")
    values_clause <- paste(rep(placeholders, nrow(batch_data)), collapse = ", ")
    query <- paste0(
      "INSERT INTO flanking_sequences (vcf_id, chromosome, position, start_position, end_position, raw_sequence, orf_nucleotide, orf_amino_acid, created_date) VALUES ",
      values_clause
    )
    
    # Flatten parameters for batch insert
    params <- list()
    for (i in 1:nrow(batch_data)) {
      params <- c(params, list(
        batch_data$vcf_id[i],
        batch_data$chromosome[i],
        batch_data$position[i],
        batch_data$start_position[i],
        batch_data$end_position[i],
        batch_data$raw_sequence[i],
        batch_data$orf_nucleotide[i],
        batch_data$orf_amino_acid[i],
        batch_data$created_date[i]
      ))
    }
    
    DBI::dbExecute(con, query, params)
  }
}

#' Export BLAST data
#' @keywords internal
.export_blast_data <- function(source_con, target_con, id_mappings, verbose) {
  
  # Check for BLAST parameters
  blast_params <- DBI::dbGetQuery(source_con, "SELECT * FROM blast_parameters ORDER BY blast_param_id")
  
  if (nrow(blast_params) == 0) {
    if (verbose) message("  - No BLAST data found in source database")
    return(id_mappings)
  }
  
  if (verbose) message("  - Exporting ", nrow(blast_params), " BLAST parameter set(s)")
  
  # Initialize BLAST ID mappings
  id_mappings$param_ids <- data.frame(old_id = integer(0), new_id = integer(0))
  id_mappings$blast_result_ids <- data.frame(old_id = integer(0), new_id = integer(0))
  
  # Export BLAST parameters
  for (i in 1:nrow(blast_params)) {
    param_record <- blast_params[i, ]
    old_param_id <- param_record$blast_param_id
    
    # Insert parameter record
    DBI::dbExecute(target_con, "
      INSERT INTO blast_parameters (blast_type, db_name, db_path, e_value, max_hits, execution_date)
      VALUES (?, ?, ?, ?, ?, ?)
    ", list(param_record$blast_type, param_record$db_name, param_record$db_path,
            param_record$e_value, param_record$max_hits, param_record$execution_date))
    
    # Get new param_id
    new_param_id <- DBI::dbGetQuery(target_con, "SELECT last_insert_rowid() as id")$id
    
    # Store parameter mapping
    id_mappings$param_ids <- rbind(id_mappings$param_ids,
                                   data.frame(old_id = old_param_id, new_id = new_param_id))
    
    # Export BLAST results for this parameter set
    blast_results <- DBI::dbGetQuery(source_con,
      "SELECT * FROM blast_results WHERE blast_param_id = ? ORDER BY blast_result_id",
      list(old_param_id))
    
    if (nrow(blast_results) > 0) {
      if (verbose) message("    - Exporting ", nrow(blast_results), " BLAST results for param set ", old_param_id)
      
      # Need to map flanking_id to new database
      blast_export <- data.frame()
      
      for (j in 1:nrow(blast_results)) {
        br_record <- blast_results[j, ]
        
        # Find corresponding flanking sequence in target database
        flanking_info <- DBI::dbGetQuery(source_con,
          "SELECT fs.chromosome, fs.position FROM flanking_sequences fs WHERE fs.flanking_id = ?",
          list(br_record$flanking_id))
        
        if (nrow(flanking_info) > 0) {
          # Find new flanking_id
          new_flanking <- DBI::dbGetQuery(target_con,
            "SELECT flanking_id FROM flanking_sequences WHERE chromosome = ? AND position = ?",
            list(flanking_info$chromosome[1], flanking_info$position[1]))
          
          if (nrow(new_flanking) > 0) {
            # Add to export data with updated IDs
            export_record <- data.frame(
              flanking_id = new_flanking$flanking_id[1],
              param_id = new_param_id,
              subject_accession = br_record$subject_accession,
              percent_identity = br_record$percent_identity,
              alignment_length = br_record$alignment_length,
              mismatches = br_record$mismatches,
              gap_opens = br_record$gap_opens,
              q_start = br_record$q_start,
              q_end = br_record$q_end,
              s_start = br_record$s_start,
              s_end = br_record$s_end,
              e_value = br_record$e_value,
              bit_score = br_record$bit_score,
              subject_title = br_record$subject_title,
              created_date = br_record$created_date
            )
            
            blast_export <- rbind(blast_export, export_record)
            
            # Store result ID mapping for annotation export
            old_result_id <- br_record$blast_result_id
            new_result_id <- nrow(blast_export)  # Will be assigned after insert
            id_mappings$blast_result_ids <- rbind(id_mappings$blast_result_ids,
                                                  data.frame(old_id = old_result_id, new_id = new_result_id))
          }
        }
      }
      
      if (nrow(blast_export) > 0) {
        # Batch insert BLAST results
        .batch_insert_blast_results(target_con, blast_export)
        
        # Update result ID mappings with actual database IDs
        result_ids <- DBI::dbGetQuery(target_con,
          "SELECT blast_result_id FROM blast_results WHERE param_id = ? ORDER BY blast_result_id",
          list(new_param_id))
        
        # Update mappings with actual IDs
        param_results <- id_mappings$blast_result_ids[id_mappings$blast_result_ids$new_id <= nrow(blast_export), ]
        for (k in 1:nrow(param_results)) {
          id_mappings$blast_result_ids$new_id[id_mappings$blast_result_ids$old_id == param_results$old_id[k]] <- result_ids$blast_result_id[k]
        }
      }
    }
  }
  
  # Export BLAST database metadata
  blast_db_metadata <- DBI::dbGetQuery(source_con, "SELECT * FROM blast_db_metadata")
  if (nrow(blast_db_metadata) > 0) {
    if (verbose) message("    - Exporting BLAST database metadata")
    
    for (i in 1:nrow(blast_db_metadata)) {
      metadata_record <- blast_db_metadata[i, ]
      DBI::dbExecute(target_con, "
        INSERT INTO blast_db_metadata (db_name, db_path, db_type, sequences, total_length, created_date)
        VALUES (?, ?, ?, ?, ?, ?)
      ", list(metadata_record$db_name, metadata_record$db_path, metadata_record$db_type,
              metadata_record$sequences, metadata_record$total_length, metadata_record$created_date))
    }
  }
  
  return(id_mappings)
}

#' Export annotation data
#' @keywords internal
.export_annotation_data <- function(source_con, target_con, id_mappings, verbose) {
  
  # Check for annotations
  annotations <- DBI::dbGetQuery(source_con, "SELECT * FROM annotations ORDER BY annotation_id")
  
  if (nrow(annotations) == 0) {
    if (verbose) message("  - No annotation data found in source database")
    return(NULL)
  }
  
  if (verbose) message("  - Exporting ", nrow(annotations), " functional annotations")
  
  # Export annotations with updated blast_result_ids
  annotation_export <- data.frame()
  
  for (i in 1:nrow(annotations)) {
    ann_record <- annotations[i, ]
    old_blast_result_id <- ann_record$blast_result_id
    
    # Find new blast_result_id from mapping
    if ("blast_result_ids" %in% names(id_mappings)) {
      new_blast_result_id <- id_mappings$blast_result_ids$new_id[id_mappings$blast_result_ids$old_id == old_blast_result_id]
      
      if (length(new_blast_result_id) > 0) {
        # Add to export data with updated blast_result_id
        export_record <- data.frame(
          blast_result_id = new_blast_result_id[1],
          uniprot_accession = ann_record$uniprot_accession,
          protein_name = ann_record$protein_name,
          gene_names = ann_record$gene_names,
          organism = ann_record$organism,
          length = ann_record$length,
          created_date = ann_record$created_date
        )
        
        annotation_export <- rbind(annotation_export, export_record)
      }
    }
  }
  
  if (nrow(annotation_export) > 0) {
    # Batch insert annotations
    .batch_insert_annotations(target_con, annotation_export)
    if (verbose) message("    - Successfully exported ", nrow(annotation_export), " annotations")
    
    # Export associated GO, KEGG, and other annotation types
    .export_functional_annotations(source_con, target_con, annotations, annotation_export, verbose)
  }
}

#' Batch insert BLAST results
#' @keywords internal
.batch_insert_blast_results <- function(con, blast_results) {
  
  # Insert in batches
  batch_size <- 500
  num_batches <- ceiling(nrow(blast_results) / batch_size)
  
  for (batch in 1:num_batches) {
    start_idx <- (batch - 1) * batch_size + 1
    end_idx <- min(batch * batch_size, nrow(blast_results))
    batch_data <- blast_results[start_idx:end_idx, ]
    
    # Create parameterized query for batch
    placeholders <- paste0("(", paste(rep("?", 14), collapse = ", "), ")")
    values_clause <- paste(rep(placeholders, nrow(batch_data)), collapse = ", ")
    query <- paste0(
      "INSERT INTO blast_results (blast_param_id, flanking_id, hit_accession, hit_description, percent_identity, alignment_length, mismatches, gap_openings, query_start, query_end, subject_start, subject_end, e_value, bit_score) VALUES ",
      values_clause
    )
    
    # Flatten parameters for batch insert
    params <- list()
    for (i in 1:nrow(batch_data)) {
      params <- c(params, list(
        batch_data$blast_param_id[i],
        batch_data$flanking_id[i],
        batch_data$hit_accession[i],
        batch_data$hit_description[i],
        batch_data$percent_identity[i],
        batch_data$alignment_length[i],
        batch_data$mismatches[i],
        batch_data$gap_openings[i],
        batch_data$query_start[i],
        batch_data$query_end[i],
        batch_data$subject_start[i],
        batch_data$subject_end[i],
        batch_data$e_value[i],
        batch_data$bit_score[i]
      ))
    }
    
    DBI::dbExecute(con, query, params)
  }
}

#' Batch insert annotations
#' @keywords internal
.batch_insert_annotations <- function(con, annotations) {
  
  # Insert in batches
  batch_size <- 500
  num_batches <- ceiling(nrow(annotations) / batch_size)
  
  for (batch in 1:num_batches) {
    start_idx <- (batch - 1) * batch_size + 1
    end_idx <- min(batch * batch_size, nrow(annotations))
    batch_data <- annotations[start_idx:end_idx, ]
    
    # Create parameterized query for batch
    placeholders <- paste0("(", paste(rep("?", 5), collapse = ", "), ")")
    values_clause <- paste(rep(placeholders, nrow(batch_data)), collapse = ", ")
    query <- paste0(
      "INSERT INTO annotations (blast_result_id, uniprot_accession, entry_name, gene_names, retrieval_date) VALUES ",
      values_clause
    )
    
    # Flatten parameters for batch insert
    params <- list()
    for (i in 1:nrow(batch_data)) {
      params <- c(params, list(
        batch_data$blast_result_id[i],
        batch_data$uniprot_accession[i],
        batch_data$entry_name[i],
        batch_data$gene_names[i],
        batch_data$retrieval_date[i]
      ))
    }
    
    DBI::dbExecute(con, query, params)
  }
}

#' Export functional annotations (GO, KEGG, etc.)
#' @keywords internal
.export_functional_annotations <- function(source_con, target_con, source_annotations, target_annotations, verbose) {
  
  # Get all annotation types to export
  annotation_tables <- c("go_terms", "kegg_references", "pfam_domains", "interpro_families", "eggnog_categories")
  
  for (table_name in annotation_tables) {
    # Check if table exists and has data
    table_exists <- tryCatch({
      count <- DBI::dbGetQuery(source_con, paste0("SELECT COUNT(*) as count FROM ", table_name))$count
      count > 0
    }, error = function(e) FALSE)
    
    if (table_exists) {
      if (verbose) message("    - Exporting ", table_name, " data")
      
      # Get functional annotation data
      func_data <- DBI::dbGetQuery(source_con, paste0("SELECT * FROM ", table_name))
      
      if (nrow(func_data) > 0) {
        # Filter to only annotations that were exported
        old_annotation_ids <- source_annotations$annotation_id
        exported_func_data <- func_data[func_data$annotation_id %in% old_annotation_ids, ]
        
        if (nrow(exported_func_data) > 0) {
          # Map old annotation_ids to new ones
          for (i in 1:nrow(exported_func_data)) {
            old_ann_id <- exported_func_data$annotation_id[i]
            # Find position in source annotations
            source_pos <- which(source_annotations$annotation_id == old_ann_id)
            if (length(source_pos) > 0) {
              # Get corresponding new annotation_id
              new_ann_id <- DBI::dbGetQuery(target_con, 
                "SELECT annotation_id FROM annotations ORDER BY annotation_id LIMIT 1 OFFSET ?",
                list(source_pos[1] - 1))$annotation_id
              exported_func_data$annotation_id[i] <- new_ann_id
            }
          }
          
          # Insert functional annotation data
          .batch_insert_functional_annotations(target_con, table_name, exported_func_data)
        }
      }
    }
  }
}

#' Batch insert functional annotations
#' @keywords internal
.batch_insert_functional_annotations <- function(con, table_name, func_data) {
  
  # Define column structures for each annotation type
  table_columns <- list(
    go_terms = c("annotation_id", "go_id", "go_term", "go_category", "go_evidence"),
    kegg_references = c("annotation_id", "kegg_id", "pathway_name"),
    pfam_domains = c("annotation_id", "pfam_id", "domain_name", "match_status"),
    interpro_families = c("annotation_id", "interpro_id", "family_name"),
    eggnog_categories = c("annotation_id", "eggnog_id", "taxonomic_scope")
  )
  
  if (!table_name %in% names(table_columns)) {
    return(NULL)
  }
  
  columns <- table_columns[[table_name]]
  num_cols <- length(columns)
  
  # Insert in batches
  batch_size <- 500
  num_batches <- ceiling(nrow(func_data) / batch_size)
  
  for (batch in 1:num_batches) {
    start_idx <- (batch - 1) * batch_size + 1
    end_idx <- min(batch * batch_size, nrow(func_data))
    batch_data <- func_data[start_idx:end_idx, ]
    
    # Create parameterized query for batch
    placeholders <- paste0("(", paste(rep("?", num_cols), collapse = ", "), ")")
    values_clause <- paste(rep(placeholders, nrow(batch_data)), collapse = ", ")
    query <- paste0(
      "INSERT INTO ", table_name, " (", paste(columns, collapse = ", "), ") VALUES ",
      values_clause
    )
    
    # Flatten parameters for batch insert
    params <- list()
    for (i in 1:nrow(batch_data)) {
      for (col in columns) {
        params <- c(params, list(batch_data[[col]][i]))
      }
    }
    
    DBI::dbExecute(con, query, params)
  }
}

#' Show export summary
#' @keywords internal
.show_export_summary <- function(target_con, verbose) {
  
  if (!verbose) return(NULL)
  
  # Get table counts
  tables <- c("input_files", "vcf_data", "reference_genomes", "reference_sequences",
              "candidate_loci", "locus_statistics", "flanking_sequences", 
              "blast_parameters", "blast_results", "annotations", "go_terms", "kegg_pathways")
  
  message("\n--- Export Summary ---")
  
  for (table in tables) {
    tryCatch({
      count <- DBI::dbGetQuery(target_con, paste0("SELECT COUNT(*) as count FROM ", table))$count
      if (count > 0) {
        message("  ", table, ": ", format(count, big.mark = ","), " records")
      }
    }, error = function(e) {
      # Table might not exist or be empty
    })
  }
}