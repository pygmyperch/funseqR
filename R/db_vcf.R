#' VCF data storage and retrieval functions for funseqR
#'
#' These functions manage VCF data within the funseqR database.
#'

#' Import a VCF file into the database
#'
#' This function reads a VCF file and stores its contents in the database.
#' It can be used as an alternative to the original read_vcf function.
#'
#' @param con A database connection object.
#' @param vcf_file A character string specifying the path to the VCF file.
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE.
#'
#' @return A list containing:
#'   \item{file_id}{The ID of the registered input file.}
#'   \item{vcf_count}{The number of VCF entries imported.}
#'
#' @importFrom vcfR read.vcfR getCHROM getPOS getID getREF getALT getQUAL getFILTER getINFO
#' @importFrom DBI dbExecute dbGetQuery
#' @importFrom progress progress_bar
#' @export
import_vcf_to_db <- function(con, vcf_file, verbose = TRUE) {
  # Register input file
  if (verbose) message("Registering VCF file...")
  file_id <- register_input_file(con, vcf_file, "vcf", verbose = verbose)

  # Check if VCF data for this file already exists
  existing <- DBI::dbGetQuery(
    con,
    "SELECT COUNT(*) AS count FROM vcf_data WHERE file_id = ?",
    params = list(file_id)
  )$count

  if (existing > 0) {
    if (verbose) message("VCF data for this file already exists in the database (", existing, " entries).")
    return(list(file_id = file_id, vcf_count = existing))
  }

  # Read VCF file
  if (verbose) message("Reading VCF file...")
  vcf_data <- vcfR::read.vcfR(vcf_file, verbose = verbose)

  # Extract data
  if (verbose) message("Extracting VCF data...")
  chrom <- vcfR::getCHROM(vcf_data)
  pos <- vcfR::getPOS(vcf_data)
  id <- vcfR::getID(vcf_data)
  ref <- vcfR::getREF(vcf_data)
  alt <- vcfR::getALT(vcf_data)
  qual <- vcfR::getQUAL(vcf_data)
  filter <- vcfR::getFILTER(vcf_data)
  info <- vcfR::getINFO(vcf_data)

  # Get FORMAT and sample data if available
  format <- NULL
  sample_data <- NULL

  if (!is.null(vcf_data@gt)) {
    format <- vcf_data@gt[, 1]  # The FORMAT column is the first column in @gt

    # Combine all sample columns into JSON
    gt_data <- vcf_data@gt[, -1, drop = FALSE]  # Exclude FORMAT column

    if (ncol(gt_data) > 0) {
      sample_data <- vapply(1:nrow(gt_data), function(i) {
        row_data <- gt_data[i, , drop = TRUE]
        names(row_data) <- colnames(gt_data)
        jsonlite::toJSON(as.list(row_data))
      }, character(1))
    }
  }

  # Prepare for import
  vcf_entries <- data.frame(
    file_id = file_id,
    chromosome = chrom,
    position = pos,
    id = id,
    ref = ref,
    alt = alt,
    qual = qual,
    filter = filter,
    info = info,
    format = format,
    sample_data = sample_data,
    stringsAsFactors = FALSE
  )

  # Start transaction
  DBI::dbExecute(con, "BEGIN TRANSACTION")

  # Set up progress bar
  if (verbose) {
    pb <- progress::progress_bar$new(
      format = "[:bar] :percent ETA: :eta",
      total = nrow(vcf_entries),
      clear = FALSE,
      width = 60
    )
  }

  # Import data in batches
  batch_size <- 1000
  num_batches <- ceiling(nrow(vcf_entries) / batch_size)

  tryCatch({
    for (i in 1:num_batches) {
      start_idx <- (i - 1) * batch_size + 1
      end_idx <- min(i * batch_size, nrow(vcf_entries))
      batch <- vcf_entries[start_idx:end_idx, ]

      # Prepare statement
      params <- list()
      for (j in 1:nrow(batch)) {
        params <- c(params, list(
          batch$file_id[j],
          batch$chromosome[j],
          batch$position[j],
          batch$id[j],
          batch$ref[j],
          batch$alt[j],
          batch$qual[j],
          batch$filter[j],
          batch$info[j],
          batch$format[j],
          batch$sample_data[j]
        ))
      }

      # Execute batch insert
      placeholders <- paste0("(", paste(rep("?", 11), collapse = ", "), ")")
      query <- paste0(
        "INSERT INTO vcf_data (file_id, chromosome, position, id, ref, alt, qual, filter, info, format, sample_data) VALUES ",
        paste(rep(placeholders, nrow(batch)), collapse = ", ")
      )

      DBI::dbExecute(con, query, params = params)

      # Update progress bar
      if (verbose) {
        pb$update(end_idx / nrow(vcf_entries))
      }
    }

    # Commit transaction
    DBI::dbExecute(con, "COMMIT")

    if (verbose) message("Imported ", nrow(vcf_entries), " VCF entries.")

    return(list(file_id = file_id, vcf_count = nrow(vcf_entries)))
  }, error = function(e) {
    # Rollback transaction on error
    DBI::dbExecute(con, "ROLLBACK")
    stop("Error importing VCF data: ", e$message)
  })

  # Update analysis report if it exists
  tryCatch({
    # Get chromosome distribution for report
    chrom_summary <- DBI::dbGetQuery(
      con,
      "SELECT chromosome, COUNT(*) as variant_count FROM vcf_data 
       WHERE file_id = ? GROUP BY chromosome ORDER BY chromosome",
      params = list(file_id)
    )
    
    n_chromosomes <- nrow(chrom_summary)
    file_size_mb <- round(file.size(vcf_file) / (1024^2), 2)
    
    vcf_message <- paste0(
      "**VCF Import Completed**\n\n",
      "- **File:** ", basename(vcf_file), " (", file_size_mb, " MB)\n",
      "- **Variants imported:** ", format(nrow(vcf_entries), big.mark = ","), "\n",
      "- **Chromosomes/Scaffolds:** ", n_chromosomes, "\n",
      "- **File ID:** ", file_id
    )
    
    update_analysis_report(
      con, project_id,
      section = "vcf_import",
      message = vcf_message,
      verbose = FALSE
    )
  }, error = function(e) {
    # Silently ignore if no report exists
  })
}

#' Get VCF data from the database
#'
#' This function retrieves VCF data from the database.
#'
#' @param con A database connection object.
#' @param file_id The ID of the input file containing the VCF data.
#' @param limit The maximum number of entries to return. Default is NULL (no limit).
#' @param offset The number of entries to skip. Default is 0.
#'
#' @return A data frame containing the VCF data.
#'
#' @importFrom DBI dbGetQuery
#' @export
get_vcf_data <- function(con, file_id, limit = NULL, offset = 0) {
  # Check if file exists
  file_info <- DBI::dbGetQuery(
    con,
    "SELECT * FROM input_files WHERE file_id = ?",
    params = list(file_id)
  )

  if (nrow(file_info) == 0) {
    stop("File with ID ", file_id, " not found.")
  }

  # Build query
  query <- "SELECT * FROM vcf_data WHERE file_id = ? ORDER BY chromosome, position"

  if (!is.null(limit)) {
    query <- paste0(query, " LIMIT ", limit, " OFFSET ", offset)
  }

  # Execute query
  DBI::dbGetQuery(con, query, params = list(file_id))
}

#' Get VCF data by chromosome and position
#'
#' This function retrieves VCF data for a specific chromosome and position range.
#'
#' @param con A database connection object.
#' @param file_id The ID of the input file containing the VCF data.
#' @param chromosome The chromosome to filter by.
#' @param start_pos The start position of the range.
#' @param end_pos The end position of the range.
#'
#' @return A data frame containing the VCF data.
#'
#' @importFrom DBI dbGetQuery
#' @export
get_vcf_by_position <- function(con, file_id, chromosome, start_pos, end_pos) {
  # Check if file exists
  file_info <- DBI::dbGetQuery(
    con,
    "SELECT * FROM input_files WHERE file_id = ?",
    params = list(file_id)
  )

  if (nrow(file_info) == 0) {
    stop("File with ID ", file_id, " not found.")
  }

  # Execute query
  DBI::dbGetQuery(
    con,
    "SELECT * FROM vcf_data
     WHERE file_id = ? AND chromosome = ? AND position >= ? AND position <= ?
     ORDER BY position",
    params = list(file_id, chromosome, start_pos, end_pos)
  )
}

#' Convert VCF data in the database to BED format
#'
#' This function generates BED format data from VCF data stored in the database.
#' It can be used as an alternative to the original vcf2bed function.
#'
#' @param con A database connection object.
#' @param file_id The ID of the input file containing the VCF data.
#' @param output_file Optional. The file path to the output BED file. If not supplied, the function will only return the data frame.
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE.
#'
#' @return A data frame with columns chromosome, start_position, and end_position.
#'
#' @importFrom DBI dbGetQuery
#' @export
vcf2bed_db <- function(con, file_id, output_file = NULL, verbose = TRUE) {
  if (verbose) message("Retrieving VCF data from database...")

  # Get VCF data
  vcf_data <- DBI::dbGetQuery(
    con,
    "SELECT chromosome, position FROM vcf_data WHERE file_id = ? ORDER BY chromosome, position",
    params = list(file_id)
  )

  if (nrow(vcf_data) == 0) {
    stop("No VCF data found for file ID ", file_id)
  }

  # Create BED data
  bed_data <- data.frame(
    chrom = vcf_data$chromosome,
    start = vcf_data$position - 1,  # BED format is 0-based
    end = vcf_data$position,
    stringsAsFactors = FALSE
  )

  # Write to file if requested
  if (!is.null(output_file)) {
    if (verbose) message("Writing BED data to file: ", output_file)
    write.table(bed_data, output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    if (verbose) message("BED file created successfully!")
  }

  return(bed_data)
}

#' Get a count of VCF entries in the database
#'
#' @param con A database connection object.
#' @param file_id The ID of the input file. If NULL, count all entries. Default is NULL.
#' @param project_id The ID of the project. If NULL, count all entries. Default is NULL.
#'
#' @return The number of VCF entries.
#'
#' @importFrom DBI dbGetQuery
#' @export
count_vcf_entries <- function(con, file_id = NULL, project_id = NULL) {
  if (!is.null(file_id)) {
    # Count entries for a specific file
    DBI::dbGetQuery(
      con,
      "SELECT COUNT(*) AS count FROM vcf_data WHERE file_id = ?",
      params = list(file_id)
    )$count
  } else if (!is.null(project_id)) {
    # Count entries for a specific project
    DBI::dbGetQuery(
      con,
      "SELECT COUNT(*) AS count FROM vcf_data v
       JOIN input_files f ON v.file_id = f.file_id
       WHERE f.project_id = ?",
      params = list(project_id)
    )$count
  } else {
    # Count all entries
    DBI::dbGetQuery(con, "SELECT COUNT(*) AS count FROM vcf_data")$count
  }
}

#' Define or discover main chromosomes for visualization consolidation
#'
#' This function serves two purposes: discovering available chromosomes in the database
#' or defining which chromosomes should be kept separate in visualizations. All other
#' chromosomes/scaffolds will be grouped into a pseudo-chromosome "US" (unmapped scaffolds).
#'
#' @param con A database connection object
#' @param main_chromosomes Character vector. Chromosome names to keep separate. 
#'   If NULL (default), function returns all available chromosome names for discovery
#' @param verbose Logical. If TRUE, print informational messages. Default is TRUE
#'
#' @return Character vector. If main_chromosomes is NULL, returns all available 
#'   chromosome names. If main_chromosomes is provided, returns them invisibly after storing
#'
#' @details
#' **Discovery mode (no parameters):** Returns all unique chromosome names found in VCF data,
#' useful for users to see what's available before choosing main chromosomes.
#' 
#' **Configuration mode (with main_chromosomes):** Stores the specified chromosome names
#' as the "main chromosomes" for this database. These will be kept separate in visualizations
#' while all others are grouped as "US".
#'
#' @examples
#' \dontrun{
#' # Discovery - see all available chromosomes
#' con <- connect_funseq_db("my_analysis.db")
#' all_chroms <- define_chromosomes(con)
#' print(all_chroms)
#' # [1] "LG1" "LG2" "LG3" "scaffold_001" "scaffold_002" ...
#' 
#' # Configuration - define main chromosomes
#' define_chromosomes(con, c("LG1", "LG2", "LG3", "LG4", "LG5"))
#' 
#' # Now visualizations will show: LG1, LG2, LG3, LG4, LG5, US
#' get_database_summary(con)  # Uses consolidated view by default
#' }
#'
#' @importFrom DBI dbGetQuery dbExecute
#' @export
define_chromosomes <- function(con, main_chromosomes = NULL, verbose = TRUE) {
  
  if (!DBI::dbIsValid(con)) {
    stop("Invalid database connection")
  }
  
  if (is.null(main_chromosomes)) {
    # Discovery mode: return all unique chromosome names
    tryCatch({
      chromosomes <- DBI::dbGetQuery(con, 
        "SELECT DISTINCT chromosome FROM vcf_data ORDER BY chromosome")$chromosome
      
      if (length(chromosomes) == 0) {
        if (verbose) message("No VCF data found in database. Import VCF files first.")
        return(character(0))
      }
      
      if (verbose) {
        cat("=== Available Chromosomes/Scaffolds ===", "\n")
        cat("Found", length(chromosomes), "chromosomes/scaffolds:\n\n")
        
        # Show in columns for better readability
        if (length(chromosomes) <= 20) {
          cat(paste(chromosomes, collapse = ", "), "\n")
        } else {
          # Show first 15, then summary
          cat(paste(chromosomes[1:15], collapse = ", "), "\n")
          cat("... and", length(chromosomes) - 15, "more\n")
        }
        
        cat("\nUsage:\n")
        cat("  define_chromosomes(con, c('LG1', 'LG2', 'LG3')) to set main chromosomes\n")
        cat("  Other scaffolds will be grouped as 'US' in visualizations\n")
      }
      
      return(chromosomes)
      
    }, error = function(e) {
      stop("Error retrieving chromosome names: ", e$message)
    })
    
  } else {
    # Configuration mode: store main chromosomes
    if (!is.character(main_chromosomes) || length(main_chromosomes) == 0) {
      stop("main_chromosomes must be a non-empty character vector")
    }
    
    # Verify these chromosomes exist
    existing_chroms <- DBI::dbGetQuery(con, 
      "SELECT DISTINCT chromosome FROM vcf_data")$chromosome
    
    missing_chroms <- setdiff(main_chromosomes, existing_chroms)
    if (length(missing_chroms) > 0) {
      warning("Some specified chromosomes not found in database: ", 
              paste(missing_chroms, collapse = ", "))
    }
    
    # Store in database metadata
    store_main_chromosomes(con, main_chromosomes)
    
    if (verbose) {
      cat("=== Main Chromosomes Defined ===", "\n")
      cat("Main chromosomes:", paste(main_chromosomes, collapse = ", "), "\n")
      cat("All other scaffolds will be grouped as 'US' in visualizations\n")
      
      valid_chroms <- intersect(main_chromosomes, existing_chroms)
      if (length(valid_chroms) > 0) {
        cat("Confirmed in database:", paste(valid_chroms, collapse = ", "), "\n")
      }
    }
    
    return(invisible(main_chromosomes))
  }
}

#' Store main chromosomes in database metadata
#' @param con Database connection
#' @param chromosomes Character vector of chromosome names
store_main_chromosomes <- function(con, chromosomes) {
  
  # Ensure metadata table exists
  tables <- DBI::dbListTables(con)
  if (!"metadata" %in% tables) {
    DBI::dbExecute(con, "
      CREATE TABLE metadata (
        key TEXT PRIMARY KEY,
        value TEXT NOT NULL
      )
    ")
  }
  
  # Store as JSON string
  chrom_json <- jsonlite::toJSON(chromosomes)
  
  # Insert or update
  existing <- DBI::dbGetQuery(con, 
    "SELECT value FROM metadata WHERE key = 'main_chromosomes'")
  
  if (nrow(existing) > 0) {
    DBI::dbExecute(con, 
      "UPDATE metadata SET value = ? WHERE key = 'main_chromosomes'",
      params = list(chrom_json))
  } else {
    DBI::dbExecute(con, 
      "INSERT INTO metadata (key, value) VALUES ('main_chromosomes', ?)",
      params = list(chrom_json))
  }
}

#' Retrieve main chromosomes from database metadata
#' @param con Database connection
#' @return Character vector of main chromosome names, or NULL if not defined
get_main_chromosomes <- function(con) {
  
  tables <- DBI::dbListTables(con)
  if (!"metadata" %in% tables) {
    return(NULL)
  }
  
  tryCatch({
    result <- DBI::dbGetQuery(con, 
      "SELECT value FROM metadata WHERE key = 'main_chromosomes'")
    
    if (nrow(result) == 0) {
      return(NULL)
    }
    
    # Parse JSON
    chromosomes <- jsonlite::fromJSON(result$value[1])
    
    # Ensure it's a character vector
    if (is.character(chromosomes)) {
      return(chromosomes)
    } else {
      return(NULL)
    }
    
  }, error = function(e) {
    return(NULL)
  })
}

#' Consolidate scaffolds into a pseudo-chromosome for visualization
#'
#' This function takes VCF data and consolidates small scaffolds into a single
#' pseudo-chromosome called "US" (unmapped scaffolds) while retaining specified
#' main chromosomes. This is useful for visualization and analysis when you have
#' many small scaffolds that clutter plots.
#'
#' @param vcf_data Data frame. VCF data with at least a 'chromosome' column
#' @param keep_chromosomes Character vector. Chromosome names to keep separate
#' @param pseudo_name Character. Name for the consolidated pseudo-chromosome. Default is "US"
#' @param verbose Logical. If TRUE, print consolidation information. Default is TRUE
#'
#' @return Data frame. VCF data with consolidated chromosome names
#'
#' @details
#' All chromosomes not specified in `keep_chromosomes` will be renamed to the
#' `pseudo_name` value. The original data is not modified - a new data frame
#' is returned with updated chromosome names.
#'
#' @examples
#' \dontrun{
#' # Get VCF data from database
#' vcf_data <- get_vcf_data(con, file_id = 1)
#' 
#' # Consolidate all except main linkage groups
#' main_chroms <- c("LG1", "LG2", "LG3", "LG4", "LG5")
#' consolidated <- consolidate_scaffolds(vcf_data, main_chroms)
#' 
#' # Check the result
#' table(consolidated$chromosome)
#' }
#'
#' @export
consolidate_scaffolds <- function(vcf_data, keep_chromosomes, pseudo_name = "US", verbose = FALSE) {
  
  if (!is.data.frame(vcf_data)) {
    stop("vcf_data must be a data frame")
  }
  
  if (!"chromosome" %in% names(vcf_data)) {
    stop("vcf_data must contain a 'chromosome' column")
  }
  
  if (!is.character(keep_chromosomes) || length(keep_chromosomes) == 0) {
    stop("keep_chromosomes must be a non-empty character vector")
  }
  
  if (nrow(vcf_data) == 0) {
    if (verbose) message("No VCF data to consolidate")
    return(vcf_data)
  }
  
  # Get original chromosome counts
  original_chroms <- table(vcf_data$chromosome)
  
  # Create copy of data
  consolidated_data <- vcf_data
  
  # Identify chromosomes to consolidate
  all_chroms <- unique(vcf_data$chromosome)
  consolidate_chroms <- setdiff(all_chroms, keep_chromosomes)
  
  if (length(consolidate_chroms) == 0) {
    if (verbose) message("No chromosomes to consolidate - all are in keep_chromosomes list")
    return(vcf_data)
  }
  
  # Replace chromosome names
  consolidated_data$chromosome[consolidated_data$chromosome %in% consolidate_chroms] <- pseudo_name
  
  if (verbose) {
    cat("=== Scaffold Consolidation Summary ===", "\n")
    cat("Kept separate:", length(keep_chromosomes), "chromosomes\n")
    cat("Consolidated:", length(consolidate_chroms), "scaffolds into", pseudo_name, "\n")
    
    # Count variants before and after
    kept_variants <- sum(original_chroms[names(original_chroms) %in% keep_chromosomes])
    consolidated_variants <- sum(original_chroms[names(original_chroms) %in% consolidate_chroms])
    
    cat("\nVariant distribution:\n")
    cat("  Kept chromosomes:", format(kept_variants, big.mark = ","), "variants\n")
    cat("  Consolidated (", pseudo_name, "):", format(consolidated_variants, big.mark = ","), "variants\n")
    
    # Show final chromosome distribution
    final_chroms <- table(consolidated_data$chromosome)
    cat("\nFinal chromosome counts:\n")
    for (i in seq_along(final_chroms)) {
      cat(sprintf("  %-10s: %s variants\n", names(final_chroms)[i], format(final_chroms[i], big.mark = ",")))
    }
  }
  
  return(consolidated_data)
}

#' Delete VCF data from the database
#'
#' @param con A database connection object.
#' @param file_id The ID of the input file containing the VCF data to delete.
#' @param confirm Logical. If TRUE, user confirmation will be required. Default is TRUE.
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE.
#'
#' @return Logical. TRUE if the deletion was successful.
#'
#' @importFrom DBI dbExecute dbGetQuery
#' @export
delete_vcf_data <- function(con, file_id, confirm = TRUE, verbose = TRUE) {
  # Check if file exists
  file_info <- DBI::dbGetQuery(
    con,
    "SELECT * FROM input_files WHERE file_id = ?",
    params = list(file_id)
  )

  if (nrow(file_info) == 0) {
    stop("File with ID ", file_id, " not found.")
  }

  # Count VCF entries
  count <- DBI::dbGetQuery(
    con,
    "SELECT COUNT(*) AS count FROM vcf_data WHERE file_id = ?",
    params = list(file_id)
  )$count

  if (count == 0) {
    if (verbose) message("No VCF data found for file ID ", file_id)
    return(TRUE)
  }

  # Confirm deletion
  if (confirm) {
    answer <- readline(paste0("Are you sure you want to delete ", count,
                             " VCF entries for file '", file_info$file_name,
                             "' (ID: ", file_id, ")? (y/n): "))

    if (tolower(answer) != "y") {
      message("VCF data deletion cancelled.")
      return(FALSE)
    }
  }

  # Start transaction
  DBI::dbExecute(con, "BEGIN TRANSACTION")

  tryCatch({
    # Delete associated flanking sequences
    if (verbose) message("Deleting associated flanking sequences...")

    # Get vcf_ids
    vcf_ids <- DBI::dbGetQuery(
      con,
      "SELECT vcf_id FROM vcf_data WHERE file_id = ?",
      params = list(file_id)
    )$vcf_id

    if (length(vcf_ids) > 0) {
      # Delete flanking sequences
      DBI::dbExecute(
        con,
        paste0("DELETE FROM flanking_sequences WHERE vcf_id IN (",
               paste(vcf_ids, collapse = ","), ")")
      )
    }

    # Delete VCF data
    if (verbose) message("Deleting VCF data...")
    deleted <- DBI::dbExecute(
      con,
      "DELETE FROM vcf_data WHERE file_id = ?",
      params = list(file_id)
    )

    # Commit transaction
    DBI::dbExecute(con, "COMMIT")

    if (verbose) message("Deleted ", deleted, " VCF entries.")

    return(TRUE)
  }, error = function(e) {
    # Rollback transaction on error
    DBI::dbExecute(con, "ROLLBACK")
    stop("Error deleting VCF data: ", e$message)
  })
}
