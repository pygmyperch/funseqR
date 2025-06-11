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
