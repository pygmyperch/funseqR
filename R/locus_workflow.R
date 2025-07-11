#' Locus Statistics and Candidate Definition Workflow
#'
#' Unified workflow for defining locus statistics and candidate loci in the streamlined
#' one-database-per-project approach.
#'

#' Define locus statistics with optional candidate threshold
#'
#' Stores statistical values for genomic loci and optionally defines candidate loci
#' based on a threshold. This unified function replaces separate statistics and
#' candidate definition steps in the streamlined workflow.
#'
#' @param con Database connection object
#' @param statistics Various formats: VCF file path, file ID, or data frame with coordinates and statistics
#' @param candidate_threshold Numeric. Optional threshold for defining candidate loci. 
#'   If provided, loci with statistics >= threshold become candidates. Default is NULL
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return List containing:
#' \itemize{
#'   \item statistics_stored: Number of locus statistics stored
#'   \item candidates_defined: Number of candidate loci defined (if threshold provided)
#'   \item import_info: Information about data import (if VCF file was imported)
#' }
#'
#' @details
#' This function implements the streamlined "one database = one analysis project" approach:
#' 
#' \strong{Input Formats:}
#' \itemize{
#'   \item \strong{VCF file path}: Imports VCF and matches to provided statistics vector
#'   \item \strong{File ID}: Uses existing VCF data and matches to statistics vector  
#'   \item \strong{Data frame}: Must contain 'chromosome', 'position', and 'statistic' columns
#' }
#' 
#' \strong{Statistics Storage:}
#' Stores locus-level statistics in the locus_statistics table with:
#' \itemize{
#'   \item Chromosome and position coordinates
#'   \item Statistical values (p-values, q-values, etc.)
#'   \item Creation timestamp
#'   \item Unique constraint prevents duplicates
#' }
#' 
#' \strong{Candidate Definition:}
#' When candidate_threshold is provided:
#' \itemize{
#'   \item Loci with statistic >= threshold become candidates
#'   \item Stored in candidate_loci table with method "threshold"
#'   \item Enables downstream enrichment analysis
#' }
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#' 
#' # Define statistics only
#' results <- define_locus_statistics(
#'   con, 
#'   statistics = data.frame(
#'     chromosome = c("LG1", "LG2", "LG1"),
#'     position = c(1000, 2000, 3000),
#'     statistic = c(0.001, 0.05, 0.02)
#'   )
#' )
#' 
#' # Define statistics and candidates in one step
#' results <- define_locus_statistics(
#'   con,
#'   statistics = my_pvalues_df,
#'   candidate_threshold = 0.01  # p-value threshold
#' )
#' 
#' # Import VCF and define statistics/candidates  
#' results <- define_locus_statistics(
#'   con,
#'   statistics = "analysis.vcf",  # Will import VCF
#'   candidate_threshold = 0.01
#' )
#' }
#'
#' @export
define_locus_statistics <- function(con, statistics, candidate_threshold = NULL, verbose = TRUE) {
  
  if (verbose) message("=== Defining Locus Statistics ===")
  
  # Initialize return object
  result <- list(
    statistics_stored = 0,
    candidates_defined = 0,
    import_info = NULL
  )
  
  # Process input statistics into standardized format
  if (verbose) message("Processing input statistics...")
  
  if (is.character(statistics) && length(statistics) == 1) {
    # Case 1: VCF file path - import and use
    if (file.exists(statistics)) {
      if (verbose) message("  - Importing VCF file: ", statistics)
      import_result <- import_vcf(con, statistics)
      result$import_info <- import_result
      
      # Get coordinates from imported VCF
      vcf_coords <- DBI::dbGetQuery(con, "
        SELECT vcf_id, chromosome, position 
        FROM vcf_data 
        WHERE file_id = ? 
        ORDER BY vcf_id
      ", list(import_result$file_id))
      
      stop("VCF file import requires corresponding statistics vector. ",
           "Please provide statistics as data frame with coordinates.")
      
    } else {
      stop("File not found: ", statistics)
    }
    
  } else if (is.numeric(statistics) && length(statistics) == 1) {
    # Case 2: File ID - use existing VCF data
    file_id <- statistics
    if (verbose) message("  - Using existing VCF file ID: ", file_id)
    
    # Get coordinates from existing VCF
    vcf_coords <- DBI::dbGetQuery(con, "
      SELECT vcf_id, chromosome, position 
      FROM vcf_data 
      WHERE file_id = ? 
      ORDER BY vcf_id
    ", list(file_id))
    
    if (nrow(vcf_coords) == 0) {
      stop("No VCF data found for file_id: ", file_id)
    }
    
    stop("File ID input requires corresponding statistics vector. ",
         "Please provide statistics as data frame with coordinates.")
    
  } else if (is.data.frame(statistics)) {
    # Case 3: Data frame with coordinates and statistics
    if (verbose) message("  - Using provided data frame")
    
    # Validate required columns
    required_cols <- c("chromosome", "position", "statistic")
    missing_cols <- required_cols[!required_cols %in% colnames(statistics)]
    if (length(missing_cols) > 0) {
      stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }
    
    # Use provided data frame
    locus_data <- statistics[, required_cols, drop = FALSE]
    
  } else {
    stop("Invalid statistics input. Provide VCF file path, file ID, or data frame with coordinates.")
  }
  
  # Validate and clean locus data
  if (verbose) message("  - Validating ", nrow(locus_data), " loci")
  
  # Remove NA values
  na_rows <- is.na(locus_data$chromosome) | is.na(locus_data$position) | is.na(locus_data$statistic)
  if (sum(na_rows) > 0) {
    if (verbose) message("  - Removing ", sum(na_rows), " rows with missing values")
    locus_data <- locus_data[!na_rows, ]
  }
  
  # Remove duplicate coordinates (keep first occurrence)
  duplicated_coords <- duplicated(paste(locus_data$chromosome, locus_data$position))
  if (sum(duplicated_coords) > 0) {
    if (verbose) message("  - Removing ", sum(duplicated_coords), " duplicate coordinates")
    locus_data <- locus_data[!duplicated_coords, ]
  }
  
  if (nrow(locus_data) == 0) {
    stop("No valid loci remaining after cleaning")
  }
  
  # Store locus statistics
  if (verbose) message("Storing locus statistics in database...")
  
  # Clear existing statistics (one database = one analysis)
  existing_count <- DBI::dbGetQuery(con, "SELECT COUNT(*) as count FROM locus_statistics")$count
  if (existing_count > 0) {
    if (verbose) message("  - Clearing ", existing_count, " existing statistics")
    DBI::dbExecute(con, "DELETE FROM locus_statistics")
  }
  
  # Prepare data for insertion
  insert_data <- data.frame(
    chromosome = locus_data$chromosome,
    position = locus_data$position,
    statistic = locus_data$statistic,
    created_date = Sys.time(),
    stringsAsFactors = FALSE
  )
  
  # Insert statistics using prepared statement for efficiency
  stmt <- DBI::dbPrepareStatement(con, "
    INSERT INTO locus_statistics (chromosome, position, statistic, created_date)
    VALUES (?, ?, ?, ?)
  ")
  
  tryCatch({
    DBI::dbWithTransaction(con, {
      for (i in 1:nrow(insert_data)) {
        DBI::dbExecuteStatement(stmt, list(
          insert_data$chromosome[i],
          insert_data$position[i], 
          insert_data$statistic[i],
          as.character(insert_data$created_date[i])
        ))
      }
    })
    result$statistics_stored <- nrow(insert_data)
    if (verbose) message("  - Stored ", result$statistics_stored, " locus statistics")
    
  }, finally = {
    DBI::dbClearResult(stmt)
  })
  
  # Optionally define candidate loci based on threshold
  if (!is.null(candidate_threshold)) {
    if (verbose) message("Defining candidate loci with threshold: ", candidate_threshold)
    
    # Clear existing candidates (one database = one analysis)  
    existing_candidates <- DBI::dbGetQuery(con, "SELECT COUNT(*) as count FROM candidate_loci")$count
    if (existing_candidates > 0) {
      if (verbose) message("  - Clearing ", existing_candidates, " existing candidates")
      DBI::dbExecute(con, "DELETE FROM candidate_loci")
    }
    
    # Identify candidates based on threshold
    candidates <- locus_data[locus_data$statistic >= candidate_threshold, ]
    
    if (nrow(candidates) > 0) {
      # Prepare candidate data for insertion
      candidate_insert <- data.frame(
        chromosome = candidates$chromosome,
        position = candidates$position,
        method = "threshold",
        threshold = candidate_threshold,
        created_date = Sys.time(),
        stringsAsFactors = FALSE
      )
      
      # Insert candidates
      candidate_stmt <- DBI::dbPrepareStatement(con, "
        INSERT INTO candidate_loci (chromosome, position, method, threshold, created_date)
        VALUES (?, ?, ?, ?, ?)
      ")
      
      tryCatch({
        DBI::dbWithTransaction(con, {
          for (i in 1:nrow(candidate_insert)) {
            DBI::dbExecuteStatement(candidate_stmt, list(
              candidate_insert$chromosome[i],
              candidate_insert$position[i],
              candidate_insert$method[i],
              candidate_insert$threshold[i],
              as.character(candidate_insert$created_date[i])
            ))
          }
        })
        result$candidates_defined <- nrow(candidate_insert)
        if (verbose) message("  - Defined ", result$candidates_defined, " candidate loci")
        
      }, finally = {
        DBI::dbClearResult(candidate_stmt)
      })
      
    } else {
      if (verbose) message("  - No loci meet threshold criteria")
    }
  }
  
  # Summary
  if (verbose) {
    message("=== Locus Statistics Definition Complete ===")
    message("Statistics stored: ", result$statistics_stored)
    if (!is.null(candidate_threshold)) {
      message("Candidates defined: ", result$candidates_defined)
      if (result$candidates_defined > 0) {
        pct_candidates <- round(100 * result$candidates_defined / result$statistics_stored, 1)
        message("Candidate percentage: ", pct_candidates, "%")
      }
    }
  }
  
  return(result)
}

#' Define candidate loci using non-threshold methods
#'
#' Defines candidate loci using methods other than statistical thresholds,
#' such as genomic regions, gene lists, or custom criteria.
#'
#' @param con Database connection object  
#' @param loci Various formats: VCF file, BED file, data frame, or file ID
#' @param method Character. Method used for candidate definition (e.g., "genomic_region", "gene_list", "custom")
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return List containing information about candidates defined
#'
#' @details
#' This function handles candidate definition methods that don't rely on 
#' statistical thresholds stored via define_locus_statistics(). Examples:
#' \itemize{
#'   \item Genomic regions of interest
#'   \item Loci near specific genes
#'   \item Custom biological criteria
#'   \item External candidate lists
#' }
#'
#' @examples
#' \dontrun{
#' # Define candidates from genomic regions
#' results <- define_candidate_loci(
#'   con,
#'   loci = "candidate_regions.bed",
#'   method = "genomic_region"
#' )
#' 
#' # Define candidates from coordinate list
#' results <- define_candidate_loci(
#'   con,
#'   loci = data.frame(
#'     chromosome = c("LG1", "LG2"),
#'     position = c(1000, 2000)
#'   ),
#'   method = "custom"
#' )
#' }
#'
#' @export
define_candidate_loci <- function(con, loci, method, verbose = TRUE) {
  
  if (verbose) message("=== Defining Candidate Loci ===")
  if (verbose) message("Method: ", method)
  
  # Initialize return object
  result <- list(
    candidates_defined = 0,
    import_info = NULL
  )
  
  # Process input loci into standardized format
  if (verbose) message("Processing input loci...")
  
  if (is.character(loci) && length(loci) == 1) {
    # Case 1: File path (VCF, BED, etc.)
    if (file.exists(loci)) {
      if (verbose) message("  - Processing file: ", loci)
      
      # Determine file type and process accordingly
      if (grepl("\\.vcf(\\.gz)?$", loci, ignore.case = TRUE)) {
        # VCF file - import and extract coordinates
        import_result <- import_vcf(con, loci)
        result$import_info <- import_result
        
        # Get coordinates from imported VCF
        loci_coords <- DBI::dbGetQuery(con, "
          SELECT DISTINCT chromosome, position 
          FROM vcf_data 
          WHERE file_id = ?
        ", list(import_result$file_id))
        
      } else if (grepl("\\.bed$", loci, ignore.case = TRUE)) {
        # BED file - parse coordinates
        stop("BED file parsing not yet implemented")
        
      } else {
        stop("Unsupported file type: ", loci)
      }
      
    } else {
      stop("File not found: ", loci)
    }
    
  } else if (is.numeric(loci) && length(loci) == 1) {
    # Case 2: File ID - use existing data
    file_id <- loci
    if (verbose) message("  - Using existing file ID: ", file_id)
    
    # Get coordinates from existing VCF
    loci_coords <- DBI::dbGetQuery(con, "
      SELECT DISTINCT chromosome, position 
      FROM vcf_data 
      WHERE file_id = ?
    ", list(file_id))
    
    if (nrow(loci_coords) == 0) {
      stop("No data found for file_id: ", file_id)
    }
    
  } else if (is.data.frame(loci)) {
    # Case 3: Data frame with coordinates
    if (verbose) message("  - Using provided data frame")
    
    # Validate required columns
    required_cols <- c("chromosome", "position")
    missing_cols <- required_cols[!required_cols %in% colnames(loci)]
    if (length(missing_cols) > 0) {
      stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }
    
    # Extract coordinates
    loci_coords <- loci[, required_cols, drop = FALSE]
    
  } else {
    stop("Invalid loci input. Provide file path, file ID, or data frame with coordinates.")
  }
  
  # Validate and clean coordinates
  if (verbose) message("  - Validating ", nrow(loci_coords), " loci")
  
  # Remove NA values
  na_rows <- is.na(loci_coords$chromosome) | is.na(loci_coords$position)
  if (sum(na_rows) > 0) {
    if (verbose) message("  - Removing ", sum(na_rows), " rows with missing values")
    loci_coords <- loci_coords[!na_rows, ]
  }
  
  # Remove duplicates
  loci_coords <- unique(loci_coords)
  
  if (nrow(loci_coords) == 0) {
    stop("No valid loci remaining after cleaning")
  }
  
  # Clear existing candidates (one database = one analysis)
  existing_candidates <- DBI::dbGetQuery(con, "SELECT COUNT(*) as count FROM candidate_loci")$count
  if (existing_candidates > 0) {
    if (verbose) message("  - Clearing ", existing_candidates, " existing candidates")
    DBI::dbExecute(con, "DELETE FROM candidate_loci")
  }
  
  # Store candidate loci
  if (verbose) message("Storing candidate loci in database...")
  
  # Prepare data for insertion
  candidate_insert <- data.frame(
    chromosome = loci_coords$chromosome,
    position = loci_coords$position,
    method = method,
    threshold = NA,  # No threshold for non-threshold methods
    created_date = Sys.time(),
    stringsAsFactors = FALSE
  )
  
  # Insert candidates
  stmt <- DBI::dbPrepareStatement(con, "
    INSERT INTO candidate_loci (chromosome, position, method, threshold, created_date)
    VALUES (?, ?, ?, ?, ?)
  ")
  
  tryCatch({
    DBI::dbWithTransaction(con, {
      for (i in 1:nrow(candidate_insert)) {
        DBI::dbExecuteStatement(stmt, list(
          candidate_insert$chromosome[i],
          candidate_insert$position[i],
          candidate_insert$method[i],
          candidate_insert$threshold[i],
          as.character(candidate_insert$created_date[i])
        ))
      }
    })
    result$candidates_defined <- nrow(candidate_insert)
    if (verbose) message("  - Stored ", result$candidates_defined, " candidate loci")
    
  }, finally = {
    DBI::dbClearResult(stmt)
  })
  
  # Summary
  if (verbose) {
    message("=== Candidate Loci Definition Complete ===")
    message("Method: ", method)
    message("Candidates defined: ", result$candidates_defined)
  }
  
  return(result)
}