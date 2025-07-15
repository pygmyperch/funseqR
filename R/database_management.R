#' Database Management Functions
#'
#' Centralized functions for database deletion, transaction management, and cleanup operations.
#' These functions provide safe, consistent interfaces for managing funseqR database contents.
#'

#' Update database schema to remove foreground_file_id from ORA analyses
#'
#' @param con Database connection object
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return Invisible NULL
#'
#' @details
#' Updates the ora_analyses table to remove the foreground_file_id column completely,
#' since the stored candidates workflow doesn't use file IDs for candidates.
#' This eliminates the need to store NULL values for unused fields.
#'
#' @export
update_ora_schema_for_stored_candidates <- function(con, verbose = TRUE) {
  if (verbose) message("Updating ORA schema to remove foreground_file_id...")
  
  # Check if ora_analyses table exists
  tables <- DBI::dbListTables(con)
  if (!"ora_analyses" %in% tables) {
    if (verbose) message("ora_analyses table doesn't exist, nothing to update")
    return(invisible(NULL))
  }
  
  # Get current schema
  schema_info <- DBI::dbGetQuery(con, "PRAGMA table_info(ora_analyses)")
  
  # Check if foreground_file_id column exists
  fg_col <- schema_info[schema_info$name == "foreground_file_id", ]
  if (nrow(fg_col) > 0) {
    if (verbose) message("Removing foreground_file_id column...")
    
    # Begin transaction
    DBI::dbExecute(con, "BEGIN TRANSACTION")
    
    tryCatch({
      # Create new table with correct schema (no foreground_file_id)
      DBI::dbExecute(con, "
        CREATE TABLE ora_analyses_new (
          analysis_id INTEGER PRIMARY KEY,
          background_file_id INTEGER NOT NULL,
          blast_param_id INTEGER,
          annotation_type TEXT NOT NULL,
          term_type TEXT NOT NULL,
          analysis_date TEXT NOT NULL,
          total_foreground_genes INTEGER,
          total_background_genes INTEGER,
          analysis_parameters TEXT,
          enrichment_method TEXT DEFAULT 'clusterprofiler',
          FOREIGN KEY (background_file_id) REFERENCES input_files (file_id),
          FOREIGN KEY (blast_param_id) REFERENCES blast_parameters (blast_param_id)
        )
      ")
      
      # Copy existing data (excluding foreground_file_id)
      existing_data <- DBI::dbGetQuery(con, "SELECT COUNT(*) as count FROM ora_analyses")
      if (existing_data$count > 0) {
        DBI::dbExecute(con, "
          INSERT INTO ora_analyses_new 
          (analysis_id, background_file_id, blast_param_id, annotation_type, term_type, 
           analysis_date, total_foreground_genes, total_background_genes, 
           analysis_parameters, enrichment_method)
          SELECT analysis_id, background_file_id, blast_param_id, annotation_type, term_type,
                 analysis_date, total_foreground_genes, total_background_genes,
                 analysis_parameters, enrichment_method
          FROM ora_analyses
        ")
        if (verbose) message("Copied ", existing_data$count, " existing records")
      }
      
      # Drop old table and rename new one
      DBI::dbExecute(con, "DROP TABLE ora_analyses")
      DBI::dbExecute(con, "ALTER TABLE ora_analyses_new RENAME TO ora_analyses")
      
      # Recreate indexes (no foreground_file_id index)
      DBI::dbExecute(con, "CREATE INDEX idx_ora_bg_file ON ora_analyses (background_file_id)")
      DBI::dbExecute(con, "CREATE INDEX idx_ora_blast_param ON ora_analyses (blast_param_id)")
      DBI::dbExecute(con, "CREATE INDEX idx_ora_annotation_type ON ora_analyses (annotation_type)")
      DBI::dbExecute(con, "CREATE INDEX idx_ora_term_type ON ora_analyses (term_type)")
      DBI::dbExecute(con, "CREATE INDEX idx_ora_method ON ora_analyses (enrichment_method)")
      
      # Commit transaction
      DBI::dbExecute(con, "COMMIT")
      
      if (verbose) message("Schema update completed successfully")
      
    }, error = function(e) {
      # Rollback on error
      DBI::dbExecute(con, "ROLLBACK")
      stop("Schema update failed: ", e$message)
    })
    
  } else {
    if (verbose) message("Schema already has correct structure (no foreground_file_id)")
  }
  
  return(invisible(NULL))
}

# ==========================
# TRANSACTION MANAGEMENT  
# ==========================

#' Begin a transaction in the database
#'
#' @param con A database connection object.
#'
#' @return Invisible NULL.
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#' begin_transaction(con)
#' # ... perform database operations ...
#' commit_transaction(con)
#' }
#'
#' @keywords internal
begin_transaction <- function(con) {
  DBI::dbExecute(con, "BEGIN TRANSACTION")
  return(invisible(NULL))
}

#' Commit a transaction in the database
#'
#' @param con A database connection object.
#'
#' @return Invisible NULL.
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#' begin_transaction(con)
#' # ... perform database operations ...
#' commit_transaction(con)
#' }
#'
#' @keywords internal
commit_transaction <- function(con) {
  DBI::dbExecute(con, "COMMIT")
  return(invisible(NULL))
}

#' Rollback a transaction in the database
#'
#' @param con A database connection object.
#'
#' @return Invisible NULL.
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#' begin_transaction(con)
#' # ... if error occurs ...
#' rollback_transaction(con)
#' }
#'
#' @keywords internal
rollback_transaction <- function(con) {
  DBI::dbExecute(con, "ROLLBACK")
  return(invisible(NULL))
}

# ==========================
# ORA RESULTS DELETION
# ==========================

#' Delete Over-Representation Analysis (ORA) results from the database
#'
#' Removes ORA analysis results and metadata from the database for specified analysis IDs.
#' This includes entries from both the ora_results and ora_analyses tables.
#'
#' @param con Database connection object
#' @param analysis_id Integer or vector of integers. Analysis ID(s) to delete
#' @param confirm Logical. Whether to prompt for user confirmation. Default is TRUE
#' @param verbose Logical. Whether to print progress messages. Default is TRUE
#'
#' @return Logical. TRUE if deletion was successful, FALSE if cancelled
#'
#' @details
#' This function safely removes ORA results from the database by:
#' \enumerate{
#'   \item Validating that the analysis ID(s) exist
#'   \item Counting associated records in both tables
#'   \item Prompting for user confirmation (if confirm = TRUE)
#'   \item Deleting records in proper order (ora_results first, then ora_analyses)
#'   \item Using transaction safety to ensure data integrity
#' }
#'
#' The function handles foreign key constraints properly and provides detailed
#' feedback on the deletion process.
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#' 
#' # Delete a single ORA analysis
#' delete_ora_results(con, analysis_id = 1)
#' 
#' # Delete multiple ORA analyses
#' delete_ora_results(con, analysis_id = c(1, 2, 3))
#' 
#' # Delete without confirmation (for scripts)
#' delete_ora_results(con, analysis_id = c(1, 2, 3), confirm = FALSE)
#' 
#' close_funseq_db(con)
#' }
#'
#' @keywords internal
delete_ora_results <- function(con, analysis_id, confirm = TRUE, verbose = TRUE) {
  
  # Validate input
  if (missing(analysis_id) || is.null(analysis_id)) {
    stop("analysis_id is required")
  }
  
  # Ensure analysis_id is numeric
  if (!is.numeric(analysis_id)) {
    stop("analysis_id must be numeric")
  }
  
  # Check if analysis ID(s) exist
  existing_analyses <- DBI::dbGetQuery(
    con,
    paste0("SELECT analysis_id, annotation_type, term_type, analysis_date, 
            total_foreground_genes, total_background_genes 
            FROM ora_analyses 
            WHERE analysis_id IN (", paste(analysis_id, collapse = ","), ")")
  )
  
  if (nrow(existing_analyses) == 0) {
    stop("No ORA analyses found with ID(s): ", paste(analysis_id, collapse = ", "))
  }
  
  # Check for missing IDs
  missing_ids <- setdiff(analysis_id, existing_analyses$analysis_id)
  if (length(missing_ids) > 0) {
    warning("Analysis ID(s) not found: ", paste(missing_ids, collapse = ", "))
  }
  
  # Count associated results
  results_count <- DBI::dbGetQuery(
    con,
    paste0("SELECT COUNT(*) AS count FROM ora_results 
            WHERE analysis_id IN (", paste(existing_analyses$analysis_id, collapse = ","), ")")
  )$count
  
  if (verbose) {
    message("Found ", nrow(existing_analyses), " ORA analysis(es) to delete:")
    for (i in 1:nrow(existing_analyses)) {
      row <- existing_analyses[i, ]
      message("  - ID ", row$analysis_id, ": ", row$annotation_type, " (", row$term_type, 
              ") - ", row$total_foreground_genes, " foreground genes")
    }
    message("  - Total associated results: ", results_count)
  }
  
  # Confirm deletion
  if (confirm) {
    if (length(analysis_id) == 1) {
      prompt_msg <- paste0("Are you sure you want to delete ORA analysis ID ", analysis_id,
                          " and its ", results_count, " associated results? (y/n): ")
    } else {
      prompt_msg <- paste0("Are you sure you want to delete ", nrow(existing_analyses), 
                          " ORA analyses and their ", results_count, " associated results? (y/n): ")
    }
    
    answer <- readline(prompt_msg)
    
    if (tolower(answer) != "y") {
      message("ORA results deletion cancelled.")
      return(FALSE)
    }
  }
  
  # Start transaction
  DBI::dbExecute(con, "BEGIN TRANSACTION")
  
  tryCatch({
    # Delete results first (due to foreign key constraints)
    if (results_count > 0) {
      deleted_results <- DBI::dbExecute(
        con,
        paste0("DELETE FROM ora_results WHERE analysis_id IN (", 
               paste(existing_analyses$analysis_id, collapse = ","), ")")
      )
      
      if (verbose) message("Deleted ", deleted_results, " ORA result records")
    }
    
    # Delete analyses
    deleted_analyses <- DBI::dbExecute(
      con,
      paste0("DELETE FROM ora_analyses WHERE analysis_id IN (", 
             paste(existing_analyses$analysis_id, collapse = ","), ")")
    )
    
    # Commit transaction
    DBI::dbExecute(con, "COMMIT")
    
    if (verbose) {
      message("Successfully deleted:")
      message("  - ", deleted_analyses, " ORA analysis record(s)")
      if (results_count > 0) {
        message("  - ", results_count, " associated result record(s)")
      }
    }
    
    return(TRUE)
    
  }, error = function(e) {
    # Rollback transaction on error
    DBI::dbExecute(con, "ROLLBACK")
    stop("Error deleting ORA results: ", e$message)
  })
}

# ==========================
# BLAST RESULTS DELETION
# ==========================

#' Delete BLAST results from the database
#'
#' Removes BLAST results from the database for a specific parameter ID.
#'
#' @param con Database connection object
#' @param blast_param_id Integer. BLAST parameter ID for which to delete results
#' @param confirm Logical. Whether to prompt for user confirmation. Default is TRUE
#' @param verbose Logical. Whether to print progress messages. Default is TRUE
#'
#' @return Logical. TRUE if deletion was successful, FALSE if cancelled
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#' delete_blast_results(con, blast_param_id = 1)
#' close_funseq_db(con)
#' }
#'
#' @keywords internal
delete_blast_results <- function(con, blast_param_id, confirm = TRUE, verbose = TRUE) {
  # Check if BLAST parameters exist
  params <- DBI::dbGetQuery(
    con,
    "SELECT * FROM blast_parameters WHERE blast_param_id = ?",
    params = list(blast_param_id)
  )

  if (nrow(params) == 0) {
    stop("BLAST parameters with ID ", blast_param_id, " not found.")
  }

  # Count results
  count <- DBI::dbGetQuery(
    con,
    "SELECT COUNT(*) AS count FROM blast_results WHERE blast_param_id = ?",
    params = list(blast_param_id)
  )$count

  if (count == 0) {
    if (verbose) message("No BLAST results found for parameter ID ", blast_param_id)
    return(TRUE)
  }

  # Confirm deletion
  if (confirm) {
    answer <- readline(paste0("Are you sure you want to delete ", count,
                             " BLAST results for parameter ID ", blast_param_id,
                             "? (y/n): "))

    if (tolower(answer) != "y") {
      message("BLAST results deletion cancelled.")
      return(FALSE)
    }
  }

  # Start transaction
  DBI::dbExecute(con, "BEGIN TRANSACTION")

  tryCatch({
    # Delete BLAST results
    deleted <- DBI::dbExecute(
      con,
      "DELETE FROM blast_results WHERE blast_param_id = ?",
      params = list(blast_param_id)
    )

    # Commit transaction
    DBI::dbExecute(con, "COMMIT")

    if (verbose) message("Deleted ", deleted, " BLAST results.")

    return(TRUE)
  }, error = function(e) {
    # Rollback transaction on error
    DBI::dbExecute(con, "ROLLBACK")
    stop("Error deleting BLAST results: ", e$message)
  })
}

# ==========================
# VCF DATA DELETION
# ==========================

#' Delete VCF data and associated sequences from the database
#'
#' Removes VCF data and associated flanking sequences from the database for a specific file ID.
#'
#' @param con Database connection object
#' @param file_id Integer. File ID for which to delete VCF data
#' @param confirm Logical. Whether to prompt for user confirmation. Default is TRUE
#' @param verbose Logical. Whether to print progress messages. Default is TRUE
#'
#' @return Logical. TRUE if deletion was successful, FALSE if cancelled
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#' delete_vcf_data(con, file_id = 1)
#' close_funseq_db(con)
#' }
#'
#' @keywords internal
delete_vcf_data <- function(con, file_id, confirm = TRUE, verbose = TRUE) {
  # Check if file exists
  file_info <- DBI::dbGetQuery(
    con,
    "SELECT * FROM input_files WHERE file_id = ? AND file_type = 'vcf'",
    params = list(file_id)
  )

  if (nrow(file_info) == 0) {
    stop("VCF file with ID ", file_id, " not found.")
  }

  # Count VCF records
  vcf_count <- DBI::dbGetQuery(
    con,
    "SELECT COUNT(*) AS count FROM vcf_data WHERE file_id = ?",
    params = list(file_id)
  )$count

  # Count flanking sequences
  flanking_count <- DBI::dbGetQuery(
    con,
    "SELECT COUNT(*) AS count FROM flanking_sequences 
     WHERE vcf_id IN (SELECT vcf_id FROM vcf_data WHERE file_id = ?)",
    params = list(file_id)
  )$count

  if (vcf_count == 0) {
    if (verbose) message("No VCF data found for file ID ", file_id)
    return(TRUE)
  }

  # Confirm deletion
  if (confirm) {
    answer <- readline(paste0("Are you sure you want to delete ", vcf_count,
                             " VCF records and ", flanking_count,
                             " flanking sequences for file '", file_info$file_name,
                             "'? (y/n): "))

    if (tolower(answer) != "y") {
      message("VCF data deletion cancelled.")
      return(FALSE)
    }
  }

  # Start transaction
  DBI::dbExecute(con, "BEGIN TRANSACTION")

  tryCatch({
    # Delete flanking sequences first (due to foreign key constraints)
    if (flanking_count > 0) {
      deleted_flanking <- DBI::dbExecute(
        con,
        "DELETE FROM flanking_sequences 
         WHERE vcf_id IN (SELECT vcf_id FROM vcf_data WHERE file_id = ?)",
        params = list(file_id)
      )
      if (verbose) message("Deleted ", deleted_flanking, " flanking sequences.")
    }

    # Delete VCF data
    deleted_vcf <- DBI::dbExecute(
      con,
      "DELETE FROM vcf_data WHERE file_id = ?",
      params = list(file_id)
    )

    # Commit transaction
    DBI::dbExecute(con, "COMMIT")

    if (verbose) message("Deleted ", deleted_vcf, " VCF records.")

    return(TRUE)
  }, error = function(e) {
    # Rollback transaction on error
    DBI::dbExecute(con, "ROLLBACK")
    stop("Error deleting VCF data: ", e$message)
  })
}

# ==========================
# FLANKING SEQUENCES DELETION
# ==========================

#' Delete flanking sequences from the database
#'
#' Removes flanking sequences from the database for a specific VCF file ID.
#'
#' @param con Database connection object
#' @param vcf_file_id Integer. VCF file ID for which to delete flanking sequences
#' @param confirm Logical. Whether to prompt for user confirmation. Default is TRUE
#' @param verbose Logical. Whether to print progress messages. Default is TRUE
#'
#' @return Logical. TRUE if deletion was successful, FALSE if cancelled
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#' delete_flanking_sequences(con, vcf_file_id = 1)
#' close_funseq_db(con)
#' }
#'
#' @keywords internal
delete_flanking_sequences <- function(con, vcf_file_id, confirm = TRUE, verbose = TRUE) {
  # Check if VCF file exists
  file_info <- DBI::dbGetQuery(
    con,
    "SELECT * FROM input_files WHERE file_id = ? AND file_type = 'vcf'",
    params = list(vcf_file_id)
  )

  if (nrow(file_info) == 0) {
    stop("VCF file with ID ", vcf_file_id, " not found.")
  }

  # Count flanking sequences
  count <- DBI::dbGetQuery(
    con,
    "SELECT COUNT(*) AS count FROM flanking_sequences 
     WHERE vcf_id IN (SELECT vcf_id FROM vcf_data WHERE file_id = ?)",
    params = list(vcf_file_id)
  )$count

  if (count == 0) {
    if (verbose) message("No flanking sequences found for VCF file ID ", vcf_file_id)
    return(TRUE)
  }

  # Confirm deletion
  if (confirm) {
    answer <- readline(paste0("Are you sure you want to delete ", count,
                             " flanking sequences for VCF file '", file_info$file_name,
                             "'? (y/n): "))

    if (tolower(answer) != "y") {
      message("Flanking sequences deletion cancelled.")
      return(FALSE)
    }
  }

  # Start transaction
  DBI::dbExecute(con, "BEGIN TRANSACTION")

  tryCatch({
    # Delete flanking sequences
    deleted <- DBI::dbExecute(
      con,
      "DELETE FROM flanking_sequences 
       WHERE vcf_id IN (SELECT vcf_id FROM vcf_data WHERE file_id = ?)",
      params = list(vcf_file_id)
    )

    # Commit transaction
    DBI::dbExecute(con, "COMMIT")

    if (verbose) message("Deleted ", deleted, " flanking sequences.")

    return(TRUE)
  }, error = function(e) {
    # Rollback transaction on error
    DBI::dbExecute(con, "ROLLBACK")
    stop("Error deleting flanking sequences: ", e$message)
  })
}

# ==========================
# ANNOTATIONS DELETION
# ==========================

#' Delete annotations from the database
#'
#' Removes annotations, GO terms, and KEGG references from the database for a specific BLAST parameter ID.
#'
#' @param con Database connection object
#' @param blast_param_id Integer. BLAST parameter ID for which to delete annotations
#' @param confirm Logical. Whether to prompt for user confirmation. Default is TRUE
#' @param verbose Logical. Whether to print progress messages. Default is TRUE
#'
#' @return Logical. TRUE if deletion was successful, FALSE if cancelled
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#' delete_annotations(con, blast_param_id = 1)
#' close_funseq_db(con)
#' }
#'
#' @keywords internal
delete_annotations <- function(con, blast_param_id, confirm = TRUE, verbose = TRUE) {
  # Check if BLAST parameters exist
  params <- DBI::dbGetQuery(
    con,
    "SELECT * FROM blast_parameters WHERE blast_param_id = ?",
    params = list(blast_param_id)
  )

  if (nrow(params) == 0) {
    stop("BLAST parameters with ID ", blast_param_id, " not found.")
  }

  # Count annotations and related data
  annotation_count <- DBI::dbGetQuery(
    con,
    "SELECT COUNT(*) AS count FROM annotations a
     JOIN blast_results br ON a.blast_result_id = br.blast_result_id
     WHERE br.blast_param_id = ?",
    params = list(blast_param_id)
  )$count

  go_count <- DBI::dbGetQuery(
    con,
    "SELECT COUNT(*) AS count FROM go_terms gt
     JOIN annotations a ON gt.annotation_id = a.annotation_id
     JOIN blast_results br ON a.blast_result_id = br.blast_result_id
     WHERE br.blast_param_id = ?",
    params = list(blast_param_id)
  )$count

  kegg_count <- DBI::dbGetQuery(
    con,
    "SELECT COUNT(*) AS count FROM kegg_references kr
     JOIN annotations a ON kr.annotation_id = a.annotation_id
     JOIN blast_results br ON a.blast_result_id = br.blast_result_id
     WHERE br.blast_param_id = ?",
    params = list(blast_param_id)
  )$count

  if (annotation_count == 0) {
    if (verbose) message("No annotations found for BLAST parameter ID ", blast_param_id)
    return(TRUE)
  }

  # Confirm deletion
  if (confirm) {
    answer <- readline(paste0("Are you sure you want to delete ", annotation_count,
                             " annotations (including ", go_count, " GO terms and ",
                             kegg_count, " KEGG references) for BLAST parameter ID ",
                             blast_param_id, "? (y/n): "))

    if (tolower(answer) != "y") {
      message("Annotations deletion cancelled.")
      return(FALSE)
    }
  }

  # Start transaction
  DBI::dbExecute(con, "BEGIN TRANSACTION")

  tryCatch({
    # Delete GO terms first (due to foreign key constraints)
    if (go_count > 0) {
      deleted_go <- DBI::dbExecute(
        con,
        "DELETE FROM go_terms WHERE annotation_id IN (
          SELECT a.annotation_id FROM annotations a
          JOIN blast_results br ON a.blast_result_id = br.blast_result_id
          WHERE br.blast_param_id = ?
        )",
        params = list(blast_param_id)
      )
      if (verbose) message("Deleted ", deleted_go, " GO terms.")
    }

    # Delete KEGG references
    if (kegg_count > 0) {
      deleted_kegg <- DBI::dbExecute(
        con,
        "DELETE FROM kegg_references WHERE annotation_id IN (
          SELECT a.annotation_id FROM annotations a
          JOIN blast_results br ON a.blast_result_id = br.blast_result_id
          WHERE br.blast_param_id = ?
        )",
        params = list(blast_param_id)
      )
      if (verbose) message("Deleted ", deleted_kegg, " KEGG references.")
    }

    # Delete annotations
    deleted_annotations <- DBI::dbExecute(
      con,
      "DELETE FROM annotations WHERE annotation_id IN (
        SELECT a.annotation_id FROM annotations a
        JOIN blast_results br ON a.blast_result_id = br.blast_result_id
        WHERE br.blast_param_id = ?
      )",
      params = list(blast_param_id)
    )

    # Commit transaction
    DBI::dbExecute(con, "COMMIT")

    if (verbose) message("Deleted ", deleted_annotations, " annotations.")

    return(TRUE)
  }, error = function(e) {
    # Rollback transaction on error
    DBI::dbExecute(con, "ROLLBACK")
    stop("Error deleting annotations: ", e$message)
  })
}

# ==========================
# UTILITY FUNCTIONS
# ==========================

#' Drop temporary table if it exists
#'
#' Safely removes a temporary table from the database if it exists.
#'
#' @param con Database connection object
#' @param table_name Character. Name of the temporary table to drop
#' @param verbose Logical. Whether to print progress messages. Default is FALSE
#'
#' @return Logical. TRUE if successful
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#' drop_temp_table(con, "temp_analysis_table")
#' }
#'
#' @keywords internal
drop_temp_table <- function(con, table_name, verbose = FALSE) {
  tryCatch({
    DBI::dbExecute(con, paste0("DROP TABLE IF EXISTS ", table_name))
    if (verbose) message("Dropped temporary table: ", table_name)
    return(TRUE)
  }, error = function(e) {
    if (verbose) warning("Could not drop table ", table_name, ": ", e$message)
    return(FALSE)
  })
}