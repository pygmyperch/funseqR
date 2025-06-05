#' Sequence data storage and retrieval functions for funseqR
#'
#' These functions manage reference genomes and flanking sequences within the funseqR database.
#'

#' Import a reference genome into the database
#'
#' This function reads a reference genome file (FASTA) and stores its contents in the database.
#'
#' @param con A database connection object.
#' @param project_id The ID of the project to associate with the reference genome.
#' @param genome_file A character string specifying the path to the FASTA file.
#' @param genome_name A character string specifying a name for the genome. Default is NULL (derived from file name).
#' @param genome_build A character string specifying the build version of the genome. Default is NULL.
#' @param store_sequences Logical. If TRUE, store the actual sequence data in the database. Default is TRUE.
#'   If FALSE, only sequence metadata will be stored, and the original file will need to be accessible.
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE.
#'
#' @return A list containing:
#'   \item{file_id}{The ID of the registered input file.}
#'   \item{genome_id}{The ID of the imported genome.}
#'   \item{sequence_count}{The number of sequences imported.}
#'
#' @importFrom Biostrings readDNAStringSet
#' @importFrom DBI dbExecute dbGetQuery
#' @importFrom progress progress_bar
#' @export
import_reference_to_db <- function(con, project_id, genome_file, genome_name = NULL, genome_build = NULL, 
                                store_sequences = TRUE, verbose = TRUE) {
  # Register input file
  if (verbose) message("Registering genome file...")
  file_id <- register_input_file(con, project_id, genome_file, "fasta", verbose = verbose)
  
  # Derive genome name from file if not provided
  if (is.null(genome_name)) {
    genome_name <- tools::file_path_sans_ext(basename(genome_file))
  }
  
  # Check if genome already exists
  existing <- DBI::dbGetQuery(
    con,
    "SELECT genome_id FROM reference_genomes WHERE file_id = ? AND genome_name = ?",
    params = list(file_id, genome_name)
  )
  
  if (nrow(existing) > 0) {
    genome_id <- existing$genome_id[1]
    
    if (verbose) message("Genome already exists in the database with ID ", genome_id)
    
    # Count sequences
    sequence_count <- DBI::dbGetQuery(
      con,
      "SELECT COUNT(*) AS count FROM reference_sequences WHERE genome_id = ?",
      params = list(genome_id)
    )$count
    
    return(list(file_id = file_id, genome_id = genome_id, sequence_count = sequence_count))
  }
  
  # Read the genome file
  if (verbose) message("Reading genome file...")
  genome_data <- Biostrings::readDNAStringSet(genome_file)
  
  # Add the genome to the database
  if (verbose) message("Adding genome to database...")
  DBI::dbExecute(
    con,
    "INSERT INTO reference_genomes (file_id, genome_name, genome_build)
     VALUES (?, ?, ?)",
    params = list(file_id, genome_name, genome_build)
  )
  
  # Get the genome ID
  genome_id <- DBI::dbGetQuery(
    con,
    "SELECT genome_id FROM reference_genomes WHERE file_id = ? AND genome_name = ?",
    params = list(file_id, genome_name)
  )$genome_id[1]
  
  # Add sequences
  if (verbose) {
    message("Adding ", length(genome_data), " sequences to database...")
    pb <- progress::progress_bar$new(
      format = "[:bar] :percent ETA: :eta",
      total = length(genome_data),
      clear = FALSE,
      width = 60
    )
  }
  
  # Start transaction
  DBI::dbExecute(con, "BEGIN TRANSACTION")
  
  # Process sequences in batches
  batch_size <- 100
  num_batches <- ceiling(length(genome_data) / batch_size)
  
  tryCatch({
    for (i in 1:num_batches) {
      start_idx <- (i - 1) * batch_size + 1
      end_idx <- min(i * batch_size, length(genome_data))
      batch <- genome_data[start_idx:end_idx]
      
      for (j in 1:length(batch)) {
        seq_idx <- start_idx + j - 1
        seq_name <- names(batch)[j]
        seq_length <- length(batch[[j]])
        
        # Some FASTA headers have additional description after the name
        seq_name <- strsplit(seq_name, "\\s+")[[1]][1]
        
        # Store the sequence if requested
        seq_data <- NULL
        if (store_sequences) {
          seq_data <- as.character(batch[[j]])
        }
        
        # Add the sequence to the database
        DBI::dbExecute(
          con,
          "INSERT INTO reference_sequences (genome_id, sequence_name, sequence_length, sequence)
           VALUES (?, ?, ?, ?)",
          params = list(genome_id, seq_name, seq_length, seq_data)
        )
        
        # Update progress bar
        if (verbose) {
          pb$update(seq_idx / length(genome_data))
        }
      }
    }
    
    # Commit transaction
    DBI::dbExecute(con, "COMMIT")
    
    if (verbose) message("Imported ", length(genome_data), " sequences.")
    
    # Update analysis report if it exists
    tryCatch({
      # Calculate genome statistics
      total_length <- sum(sapply(genome_data, length))
      file_size_mb <- round(file.size(genome_file) / (1024^2), 2)
      avg_seq_length <- round(total_length / length(genome_data))
      
      genome_message <- paste0(
        "**Reference Genome Import Completed**\n\n",
        "- **File:** ", basename(genome_file), " (", file_size_mb, " MB)\n",
        "- **Genome:** ", genome_name %||% "Unknown", 
        if (!is.null(genome_build)) paste0(" (build ", genome_build, ")") else "", "\n",
        "- **Sequences imported:** ", format(length(genome_data), big.mark = ","), "\n",
        "- **Total length:** ", format(total_length, big.mark = ","), " bp\n",
        "- **Average sequence length:** ", format(avg_seq_length, big.mark = ","), " bp\n",
        "- **Genome ID:** ", genome_id
      )
      
      update_analysis_report(
        con, project_id,
        section = "genome_import",
        message = genome_message,
        verbose = FALSE
      )
    }, error = function(e) {
      # Silently ignore if no report exists
    })
    
    return(list(file_id = file_id, genome_id = genome_id, sequence_count = length(genome_data)))
  }, error = function(e) {
    # Rollback transaction on error
    DBI::dbExecute(con, "ROLLBACK")
    stop("Error importing genome: ", e$message)
  })
}

#' Get a reference genome from the database
#'
#' This function retrieves a reference genome from the database.
#'
#' @param con A database connection object.
#' @param genome_id The ID of the genome to retrieve.
#' @param as_dna_string_set Logical. If TRUE, return a DNAStringSet object. Default is TRUE.
#'   If FALSE, return a data frame with sequence data.
#' @param include_sequences Logical. If TRUE, include the actual sequence data. Default is TRUE.
#'   If FALSE, only include metadata.
#'
#' @return A DNAStringSet object or a data frame, depending on as_dna_string_set.
#'
#' @importFrom DBI dbGetQuery
#' @importFrom Biostrings DNAStringSet
#' @export
get_reference_genome <- function(con, genome_id, as_dna_string_set = TRUE, include_sequences = TRUE) {
  # Check if genome exists
  genome_info <- DBI::dbGetQuery(
    con,
    "SELECT * FROM reference_genomes WHERE genome_id = ?",
    params = list(genome_id)
  )
  
  if (nrow(genome_info) == 0) {
    stop("Genome with ID ", genome_id, " not found.")
  }
  
  # Build query
  if (include_sequences) {
    query <- "SELECT * FROM reference_sequences WHERE genome_id = ? ORDER BY sequence_id"
  } else {
    query <- "SELECT sequence_id, genome_id, sequence_name, sequence_length FROM reference_sequences WHERE genome_id = ? ORDER BY sequence_id"
  }
  
  # Get sequence data
  sequences <- DBI::dbGetQuery(con, query, params = list(genome_id))
  
  if (nrow(sequences) == 0) {
    stop("No sequences found for genome with ID ", genome_id)
  }
  
  # Check if sequences are stored in the database
  if (include_sequences && is.null(sequences$sequence[1])) {
    # Get file path and read from file
    file_info <- DBI::dbGetQuery(
      con,
      "SELECT f.file_path FROM input_files f
       JOIN reference_genomes g ON f.file_id = g.file_id
       WHERE g.genome_id = ?",
      params = list(genome_id)
    )
    
    if (nrow(file_info) == 0 || !file.exists(file_info$file_path[1])) {
      stop("Sequences not stored in database and original file not found.")
    }
    
    # Read sequences from file
    genome_data <- Biostrings::readDNAStringSet(file_info$file_path[1])
    
    # Return as requested format
    if (as_dna_string_set) {
      return(genome_data)
    } else {
      return(data.frame(
        sequence_id = 1:length(genome_data),
        genome_id = genome_id,
        sequence_name = names(genome_data),
        sequence_length = width(genome_data),
        sequence = as.character(genome_data),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Return as requested format
  if (as_dna_string_set && include_sequences) {
    # Create DNAStringSet object
    genome_data <- Biostrings::DNAStringSet(sequences$sequence)
    names(genome_data) <- sequences$sequence_name
    return(genome_data)
  } else {
    # Return as data frame
    return(sequences)
  }
}

#' Import flanking sequences into the database
#'
#' This function extracts flanking sequences for VCF data in the database and stores them.
#' It can be used as an alternative to the original extract_flanking_sequences function.
#'
#' @param con A database connection object.
#' @param vcf_file_id The ID of the input file containing the VCF data.
#' @param genome_id The ID of the reference genome to use.
#' @param flank_size The size of the flanking region to extract on each side of the SNP. Default is 10000.
#' @param chromosome Optional. Limit extraction to a specific chromosome. Default is NULL (all chromosomes).
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE.
#'
#' @return A list containing:
#'   \item{vcf_count}{The number of VCF entries processed.}
#'   \item{flanking_count}{The number of flanking sequences extracted.}
#'
#' @importFrom DBI dbExecute dbGetQuery
#' @importFrom progress progress_bar
#' @importFrom Biostrings subseq
#' @export
import_flanking_seqs_to_db <- function(con, vcf_file_id, genome_id, flank_size = 10000, 
                                 chromosome = NULL, verbose = TRUE) {
  # Check if VCF file exists
  vcf_file_info <- DBI::dbGetQuery(
    con,
    "SELECT * FROM input_files WHERE file_id = ?",
    params = list(vcf_file_id)
  )
  
  if (nrow(vcf_file_info) == 0) {
    stop("VCF file with ID ", vcf_file_id, " not found.")
  }
  
  # Check if genome exists
  genome_info <- DBI::dbGetQuery(
    con,
    "SELECT * FROM reference_genomes WHERE genome_id = ?",
    params = list(genome_id)
  )
  
  if (nrow(genome_info) == 0) {
    stop("Genome with ID ", genome_id, " not found.")
  }
  
  # Get VCF data
  if (verbose) message("Retrieving VCF data...")
  
  if (is.null(chromosome)) {
    vcf_query <- "SELECT vcf_id, chromosome, position FROM vcf_data WHERE file_id = ? ORDER BY chromosome, position"
    vcf_params <- list(vcf_file_id)
  } else {
    vcf_query <- "SELECT vcf_id, chromosome, position FROM vcf_data WHERE file_id = ? AND chromosome = ? ORDER BY position"
    vcf_params <- list(vcf_file_id, chromosome)
  }
  
  vcf_data <- DBI::dbGetQuery(con, vcf_query, params = vcf_params)
  
  if (nrow(vcf_data) == 0) {
    stop("No VCF data found for file ID ", vcf_file_id)
  }
  
  # Get reference sequences
  if (verbose) message("Retrieving reference sequences...")
  ref_sequences <- DBI::dbGetQuery(
    con,
    "SELECT sequence_id, sequence_name, sequence_length, sequence FROM reference_sequences WHERE genome_id = ?",
    params = list(genome_id)
  )
  
  if (nrow(ref_sequences) == 0) {
    stop("No reference sequences found for genome ID ", genome_id)
  }
  
  # Check if sequences are stored in the database
  if (is.null(ref_sequences$sequence[1])) {
    if (verbose) message("Sequences not stored in database, reading from file...")
    
    # Get file path
    file_path <- DBI::dbGetQuery(
      con,
      "SELECT file_path FROM input_files WHERE file_id = ?",
      params = list(genome_info$file_id)
    )$file_path[1]
    
    if (!file.exists(file_path)) {
      stop("Reference genome file not found: ", file_path)
    }
    
    # Read genome data
    genome_data <- Biostrings::readDNAStringSet(file_path)
    
    # Create a mapping from sequence names to sequences
    seq_names <- names(genome_data)
    seq_names <- sapply(strsplit(seq_names, "\\s+"), `[`, 1)  # Extract first word from FASTA headers
    
    ref_seq_map <- setNames(as.list(genome_data), seq_names)
  } else {
    # Create a mapping from sequence names to sequences
    ref_seq_map <- setNames(
      lapply(1:nrow(ref_sequences), function(i) ref_sequences$sequence[i]), 
      ref_sequences$sequence_name
    )
  }
  
  # Create a mapping from sequence names to sequence IDs
  seq_id_map <- setNames(ref_sequences$sequence_id, ref_sequences$sequence_name)
  
  # Set up progress bar
  if (verbose) {
    pb <- progress::progress_bar$new(
      format = "[:bar] :percent ETA: :eta",
      total = nrow(vcf_data),
      clear = FALSE,
      width = 60
    )
  }
  
  # Start transaction
  DBI::dbExecute(con, "BEGIN TRANSACTION")
  
  # Process VCF entries
  flanking_count <- 0
  
  tryCatch({
    for (i in 1:nrow(vcf_data)) {
      vcf_id <- vcf_data$vcf_id[i]
      chrom <- vcf_data$chromosome[i]
      pos <- vcf_data$position[i]
      
      # Check if flanking sequence already exists
      existing <- DBI::dbGetQuery(
        con,
        "SELECT flanking_id FROM flanking_sequences WHERE vcf_id = ?",
        params = list(vcf_id)
      )
      
      if (nrow(existing) > 0) {
        if (verbose) {
          pb$update(i / nrow(vcf_data))
        }
        next
      }
      
      # Check if chromosome exists in reference sequences
      if (!chrom %in% names(ref_seq_map)) {
        if (verbose) message("Chromosome ", chrom, " not found in reference genome. Skipping.")
        if (verbose) {
          pb$update(i / nrow(vcf_data))
        }
        next
      }
      
      # Get sequence ID
      seq_id <- seq_id_map[chrom]
      
      # Get reference sequence
      ref_seq <- ref_seq_map[[chrom]]
      seq_length <- nchar(ref_seq)
      
      # Extract flanking sequence
      start_pos <- max(1, pos - flank_size)
      end_pos <- min(seq_length, pos + flank_size)
      
      flanking_seq <- substr(ref_seq, start_pos, end_pos)
      
      # Store flanking sequence in database
      DBI::dbExecute(
        con,
        "INSERT INTO flanking_sequences (vcf_id, sequence_id, flank_size, start_position, end_position, sequence)
         VALUES (?, ?, ?, ?, ?, ?)",
        params = list(vcf_id, seq_id, flank_size, start_pos, end_pos, flanking_seq)
      )
      
      flanking_count <- flanking_count + 1
      
      # Update progress bar
      if (verbose) {
        pb$update(i / nrow(vcf_data))
      }
    }
    
    # Commit transaction
    DBI::dbExecute(con, "COMMIT")
    
    if (verbose) message("Extracted ", flanking_count, " flanking sequences.")
    
    # Update analysis report if it exists
    tryCatch({
      # Get project ID from VCF file
      project_info <- DBI::dbGetQuery(
        con,
        "SELECT if.project_id FROM input_files if 
         JOIN vcf_data v ON if.file_id = v.file_id 
         WHERE v.file_id = ? LIMIT 1",
        params = list(vcf_file_id)
      )
      
      if (nrow(project_info) > 0) {
        # Get genome info
        genome_info <- DBI::dbGetQuery(
          con,
          "SELECT g.genome_name, g.genome_build FROM genomes g WHERE g.genome_id = ?",
          params = list(genome_id)
        )
        
        success_rate <- round((flanking_count / nrow(vcf_data)) * 100, 1)
        
        flanking_message <- paste0(
          "**Flanking Sequence Extraction Completed**\n\n",
          "- **Reference genome:** ", genome_info$genome_name %||% "Unknown", 
          if (!is.null(genome_info$genome_build)) paste0(" (", genome_info$genome_build, ")") else "", "\n",
          "- **Flank size:** ±", format(flank_size, big.mark = ","), " bp\n",
          "- **Total variants:** ", format(nrow(vcf_data), big.mark = ","), "\n",
          "- **Successfully extracted:** ", format(flanking_count, big.mark = ","), " (", success_rate, "%)\n",
          "- **Skipped variants:** ", format(nrow(vcf_data) - flanking_count, big.mark = ",")
        )
        
        update_analysis_report(
          con, project_info$project_id[1],
          section = "flanking_extraction",
          message = flanking_message,
          verbose = FALSE
        )
      }
    }, error = function(e) {
      # Silently ignore if no report exists
    })
    
    return(list(vcf_count = nrow(vcf_data), flanking_count = flanking_count))
  }, error = function(e) {
    # Rollback transaction on error
    DBI::dbExecute(con, "ROLLBACK")
    stop("Error extracting flanking sequences: ", e$message)
  })
}

#' Get flanking sequences from the database
#'
#' @param con A database connection object.
#' @param vcf_file_id The ID of the input file containing the VCF data.
#' @param as_dna_string_set Logical. If TRUE, return a DNAStringSet object. Default is TRUE.
#'   If FALSE, return a data frame with sequence data.
#' @param chromosome Optional. Limit retrieval to a specific chromosome. Default is NULL (all chromosomes).
#'
#' @return A DNAStringSet object or a data frame, depending on as_dna_string_set.
#'
#' @importFrom DBI dbGetQuery
#' @importFrom Biostrings DNAStringSet
#' @export
get_flanking_sequences <- function(con, vcf_file_id, as_dna_string_set = TRUE, chromosome = NULL) {
  # Build query
  if (is.null(chromosome)) {
    query <- "
      SELECT f.flanking_id, v.chromosome, v.position, f.flank_size, f.start_position, f.end_position, f.sequence
      FROM flanking_sequences f
      JOIN vcf_data v ON f.vcf_id = v.vcf_id
      WHERE v.file_id = ?
      ORDER BY v.chromosome, v.position
    "
    params <- list(vcf_file_id)
  } else {
    query <- "
      SELECT f.flanking_id, v.chromosome, v.position, f.flank_size, f.start_position, f.end_position, f.sequence
      FROM flanking_sequences f
      JOIN vcf_data v ON f.vcf_id = v.vcf_id
      WHERE v.file_id = ? AND v.chromosome = ?
      ORDER BY v.position
    "
    params <- list(vcf_file_id, chromosome)
  }
  
  # Execute query
  flanking_data <- DBI::dbGetQuery(con, query, params = params)
  
  if (nrow(flanking_data) == 0) {
    stop("No flanking sequences found for VCF file ID ", vcf_file_id)
  }
  
  # Return as requested format
  if (as_dna_string_set) {
    # Create names for the sequences
    names <- paste0(
      flanking_data$chromosome, ":",
      flanking_data$position, " (",
      flanking_data$start_position, "-",
      flanking_data$end_position, ")"
    )
    
    # Create DNAStringSet object
    flanking_seqs <- Biostrings::DNAStringSet(flanking_data$sequence)
    names(flanking_seqs) <- names
    
    return(flanking_seqs)
  } else {
    return(flanking_data)
  }
}

#' Delete flanking sequences from the database
#'
#' @param con A database connection object.
#' @param vcf_file_id The ID of the input file containing the VCF data.
#' @param confirm Logical. If TRUE, user confirmation will be required. Default is TRUE.
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE.
#'
#' @return Logical. TRUE if the deletion was successful.
#'
#' @importFrom DBI dbExecute dbGetQuery
#' @export
delete_flanking_sequences <- function(con, vcf_file_id, confirm = TRUE, verbose = TRUE) {
  # Count flanking sequences
  count_query <- "
    SELECT COUNT(*) AS count
    FROM flanking_sequences f
    JOIN vcf_data v ON f.vcf_id = v.vcf_id
    WHERE v.file_id = ?
  "
  
  count <- DBI::dbGetQuery(con, count_query, params = list(vcf_file_id))$count
  
  if (count == 0) {
    if (verbose) message("No flanking sequences found for VCF file ID ", vcf_file_id)
    return(TRUE)
  }
  
  # Confirm deletion
  if (confirm) {
    answer <- readline(paste0("Are you sure you want to delete ", count, 
                             " flanking sequences for VCF file ID ", vcf_file_id, 
                             "? (y/n): "))
    
    if (tolower(answer) != "y") {
      message("Flanking sequences deletion cancelled.")
      return(FALSE)
    }
  }
  
  # Start transaction
  DBI::dbExecute(con, "BEGIN TRANSACTION")
  
  tryCatch({
    # Get vcf_ids
    vcf_ids_query <- "SELECT vcf_id FROM vcf_data WHERE file_id = ?"
    vcf_ids <- DBI::dbGetQuery(con, vcf_ids_query, params = list(vcf_file_id))$vcf_id
    
    if (length(vcf_ids) > 0) {
      # Delete flanking sequences
      delete_query <- paste0(
        "DELETE FROM flanking_sequences WHERE vcf_id IN (",
        paste(vcf_ids, collapse = ","),
        ")"
      )
      
      deleted <- DBI::dbExecute(con, delete_query)
      
      # Commit transaction
      DBI::dbExecute(con, "COMMIT")
      
      if (verbose) message("Deleted ", deleted, " flanking sequences.")
      
      return(TRUE)
    } else {
      DBI::dbExecute(con, "COMMIT")
      if (verbose) message("No VCF entries found for file ID ", vcf_file_id)
      return(TRUE)
    }
  }, error = function(e) {
    # Rollback transaction on error
    DBI::dbExecute(con, "ROLLBACK")
    stop("Error deleting flanking sequences: ", e$message)
  })
}
