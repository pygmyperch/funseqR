#' Project management functions for funseqR
#'
#' These functions manage projects within the funseqR database.
#'

#' Create a new project in the database
#'
#' This function creates a new project in the funseqR database.
#'
#' @param con A database connection object.
#' @param project_name Character string specifying the name of the project.
#' @param description Character string describing the project. Default is NULL.
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE.
#'
#' @return The ID of the newly created project.
#'
#' @importFrom DBI dbExecute dbGetQuery
#' @export
create_project <- function(con, project_name, description = NULL,
                           create_report = TRUE, report_format = "Rmd",
                           report_template = "standard", verbose = TRUE) {
  # Check if project name already exists
  existing <- DBI::dbGetQuery(
    con,
    "SELECT project_id FROM projects WHERE project_name = ?",
    params = list(project_name)
  )

  if (nrow(existing) > 0) {
    stop("A project with the name '", project_name, "' already exists.")
  }

  # Create project
  current_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  DBI::dbExecute(
    con,
    "INSERT INTO projects (project_name, description, creation_date, last_modified)
     VALUES (?, ?, ?, ?)",
    params = list(project_name, description, current_time, current_time)
  )

  # Get the ID of the newly created project
  project_id <- DBI::dbGetQuery(
    con,
    "SELECT project_id FROM projects WHERE project_name = ?",
    params = list(project_name)
  )$project_id[1]

  if (verbose) message("Created project '", project_name, "' with ID ", project_id)

  if (create_report) {
    tryCatch({
      report_path <- create_analysis_report(
        con, project_id,
        format = report_format,
        template = report_template,
        open_report = interactive(),
        verbose = verbose
      )

      if (verbose) message("Created analysis report: ", report_path)
    }, error = function(e) {
      if (verbose) message("Could not create analysis report: ", e$message)
    })
  }

  return(project_id)
}

#' Get a list of projects in the database
#'
#' @param con A database connection object.
#'
#' @return A data frame containing information about all projects in the database.
#'
#' @importFrom DBI dbGetQuery
#' @export
list_projects <- function(con) {
  DBI::dbGetQuery(con, "SELECT * FROM projects ORDER BY project_id")
}

#' Get information about a specific project
#'
#' @param con A database connection object.
#' @param project_id The ID of the project.
#'
#' @return A data frame containing information about the project.
#'
#' @importFrom DBI dbGetQuery
#' @export
get_project <- function(con, project_id) {
  project <- DBI::dbGetQuery(
    con,
    "SELECT * FROM projects WHERE project_id = ?",
    params = list(project_id)
  )

  if (nrow(project) == 0) {
    stop("Project with ID ", project_id, " not found.")
  }

  return(project)
}

#' Get the ID of a project by name
#'
#' @param con A database connection object.
#' @param project_name The name of the project.
#'
#' @return The ID of the project, or NULL if no project with the given name exists.
#'
#' @importFrom DBI dbGetQuery
#' @export
get_project_id <- function(con, project_name) {
  project <- DBI::dbGetQuery(
    con,
    "SELECT project_id FROM projects WHERE project_name = ?",
    params = list(project_name)
  )

  if (nrow(project) == 0) {
    return(NULL)
  }

  return(project$project_id[1])
}

#' Update a project in the database
#'
#' @param con A database connection object.
#' @param project_id The ID of the project to update.
#' @param project_name New name for the project. Default is NULL (no change).
#' @param description New description for the project. Default is NULL (no change).
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE.
#'
#' @return Logical. TRUE if the update was successful.
#'
#' @importFrom DBI dbExecute dbGetQuery
#' @export
update_project <- function(con, project_id, project_name = NULL, description = NULL, verbose = TRUE) {
  # Check if project exists
  project <- DBI::dbGetQuery(
    con,
    "SELECT * FROM projects WHERE project_id = ?",
    params = list(project_id)
  )

  if (nrow(project) == 0) {
    stop("Project with ID ", project_id, " not found.")
  }

  # Prepare update
  update_parts <- c()
  params <- list()

  if (!is.null(project_name)) {
    # Check if the new name conflicts with existing projects
    if (project_name != project$project_name) {
      existing <- DBI::dbGetQuery(
        con,
        "SELECT project_id FROM projects WHERE project_name = ? AND project_id != ?",
        params = list(project_name, project_id)
      )

      if (nrow(existing) > 0) {
        stop("A project with the name '", project_name, "' already exists.")
      }
    }

    update_parts <- c(update_parts, "project_name = ?")
    params <- c(params, list(project_name))
  }

  if (!is.null(description)) {
    update_parts <- c(update_parts, "description = ?")
    params <- c(params, list(description))
  }

  # Add last modified timestamp
  update_parts <- c(update_parts, "last_modified = ?")
  params <- c(params, list(format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

  # Add project_id to params
  params <- c(params, list(project_id))

  # Execute update
  if (length(update_parts) > 0) {
    query <- paste("UPDATE projects SET", paste(update_parts, collapse = ", "), "WHERE project_id = ?")
    DBI::dbExecute(con, query, params = params)

    if (verbose) message("Updated project with ID ", project_id)
  } else {
    if (verbose) message("No changes to update for project with ID ", project_id)
  }

  return(TRUE)
}

#' Delete a project from the database
#'
#' This function deletes a project and all associated data from the database.
#'
#' @param con A database connection object.
#' @param project_id The ID of the project to delete.
#' @param confirm Logical. If TRUE, user confirmation will be required. Default is TRUE.
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE.
#'
#' @return Logical. TRUE if the deletion was successful.
#'
#' @importFrom DBI dbExecute dbGetQuery
#' @export
delete_project <- function(con, project_id, confirm = TRUE, verbose = TRUE) {
  # Check if project exists
  project <- DBI::dbGetQuery(
    con,
    "SELECT project_name FROM projects WHERE project_id = ?",
    params = list(project_id)
  )

  if (nrow(project) == 0) {
    stop("Project with ID ", project_id, " not found.")
  }

  # Confirm deletion
  if (confirm) {
    answer <- readline(paste0("Are you sure you want to delete the project '",
                             project$project_name, "' (ID: ", project_id,
                             ") and all associated data? (y/n): "))

    if (tolower(answer) != "y") {
      message("Project deletion cancelled.")
      return(FALSE)
    }
  }

  # Start transaction
  DBI::dbExecute(con, "BEGIN TRANSACTION")

  tryCatch({
    # Delete all associated data
    # The order matters due to foreign key constraints

    # Delete annotations and related data
    if (verbose) message("Deleting annotations data...")

    # Get all annotation IDs
    annotation_ids <- DBI::dbGetQuery(
      con,
      "SELECT a.annotation_id FROM annotations a
       JOIN blast_results br ON a.blast_result_id = br.blast_result_id
       JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
       WHERE bp.project_id = ?",
      params = list(project_id)
    )$annotation_id

    if (length(annotation_ids) > 0) {
      # Delete GO terms
      DBI::dbExecute(
        con,
        paste0("DELETE FROM go_terms WHERE annotation_id IN (",
               paste(annotation_ids, collapse = ","), ")")
      )

      # Delete KEGG references
      DBI::dbExecute(
        con,
        paste0("DELETE FROM kegg_references WHERE annotation_id IN (",
               paste(annotation_ids, collapse = ","), ")")
      )

      # Delete annotations
      DBI::dbExecute(
        con,
        paste0("DELETE FROM annotations WHERE annotation_id IN (",
               paste(annotation_ids, collapse = ","), ")")
      )
    }

    # Delete BLAST results and parameters
    if (verbose) message("Deleting BLAST data...")

    # Get all blast parameter IDs
    blast_param_ids <- DBI::dbGetQuery(
      con,
      "SELECT blast_param_id FROM blast_parameters WHERE project_id = ?",
      params = list(project_id)
    )$blast_param_id

    if (length(blast_param_ids) > 0) {
      # Delete BLAST results
      DBI::dbExecute(
        con,
        paste0("DELETE FROM blast_results WHERE blast_param_id IN (",
               paste(blast_param_ids, collapse = ","), ")")
      )

      # Delete BLAST parameters
      DBI::dbExecute(
        con,
        paste0("DELETE FROM blast_parameters WHERE blast_param_id IN (",
               paste(blast_param_ids, collapse = ","), ")")
      )
    }

    # Delete flanking sequences
    if (verbose) message("Deleting sequence data...")

    # Get all file IDs for this project
    file_ids <- DBI::dbGetQuery(
      con,
      "SELECT file_id FROM input_files WHERE project_id = ?",
      params = list(project_id)
    )$file_id

    if (length(file_ids) > 0) {
      # Get all VCF IDs
      vcf_ids <- DBI::dbGetQuery(
        con,
        paste0("SELECT vcf_id FROM vcf_data WHERE file_id IN (",
               paste(file_ids, collapse = ","), ")")
      )$vcf_id

      if (length(vcf_ids) > 0) {
        # Delete flanking sequences
        DBI::dbExecute(
          con,
          paste0("DELETE FROM flanking_sequences WHERE vcf_id IN (",
                 paste(vcf_ids, collapse = ","), ")")
        )
      }

      # Get all genome IDs
      genome_ids <- DBI::dbGetQuery(
        con,
        paste0("SELECT genome_id FROM reference_genomes WHERE file_id IN (",
               paste(file_ids, collapse = ","), ")")
      )$genome_id

      if (length(genome_ids) > 0) {
        # Delete reference sequences
        DBI::dbExecute(
          con,
          paste0("DELETE FROM reference_sequences WHERE genome_id IN (",
                 paste(genome_ids, collapse = ","), ")")
        )

        # Delete reference genomes
        DBI::dbExecute(
          con,
          paste0("DELETE FROM reference_genomes WHERE genome_id IN (",
                 paste(genome_ids, collapse = ","), ")")
        )
      }

      # Delete VCF data
      if (length(vcf_ids) > 0) {
        DBI::dbExecute(
          con,
          paste0("DELETE FROM vcf_data WHERE file_id IN (",
                 paste(file_ids, collapse = ","), ")")
        )
      }

      # Delete input files
      DBI::dbExecute(
        con,
        paste0("DELETE FROM input_files WHERE file_id IN (",
               paste(file_ids, collapse = ","), ")")
      )
    }

    # Delete project
    if (verbose) message("Deleting project...")
    DBI::dbExecute(
      con,
      "DELETE FROM projects WHERE project_id = ?",
      params = list(project_id)
    )

    # Commit transaction
    DBI::dbExecute(con, "COMMIT")

    if (verbose) message("Project '", project$project_name, "' (ID: ", project_id, ") deleted successfully.")

    return(TRUE)
  }, error = function(e) {
    # Rollback transaction on error
    DBI::dbExecute(con, "ROLLBACK")
    stop("Error deleting project: ", e$message)
  })
}

#' Register an input file in the database
#'
#' This function adds information about an input file to the database.
#'
#' @param con A database connection object.
#' @param project_id The ID of the project to associate with the file.
#' @param file_path Character string specifying the path to the file.
#' @param file_type Character string specifying the type of the file (e.g., "vcf", "fasta").
#' @param calculate_hash Logical. If TRUE, calculate and store a hash of the file. Default is TRUE.
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE.
#'
#' @return The ID of the newly registered file.
#'
#' @importFrom DBI dbExecute dbGetQuery
#' @export
register_input_file <- function(con, project_id, file_path, file_type, calculate_hash = TRUE, verbose = TRUE) {
  # Check if project exists
  project <- DBI::dbGetQuery(
    con,
    "SELECT project_id FROM projects WHERE project_id = ?",
    params = list(project_id)
  )

  if (nrow(project) == 0) {
    stop("Project with ID ", project_id, " not found.")
  }

  # Check if file exists
  if (!file.exists(file_path)) {
    stop("File does not exist: ", file_path)
  }

  # Get file name
  file_name <- basename(file_path)

  # Calculate hash if requested
  file_hash <- NULL
  if (calculate_hash) {
    if (verbose) message("Calculating file hash...")
    file_hash <- digest::digest(file_path, algo = "sha256", file = TRUE)
  }

  # Check if this file has already been registered for this project
  existing <- DBI::dbGetQuery(
    con,
    "SELECT file_id FROM input_files
     WHERE project_id = ? AND file_name = ? AND file_path = ?",
    params = list(project_id, file_name, file_path)
  )

  if (nrow(existing) > 0) {
    if (verbose) message("File already registered with ID ", existing$file_id[1])
    return(existing$file_id[1])
  }

  # Register file
  current_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  DBI::dbExecute(
    con,
    "INSERT INTO input_files (project_id, file_type, file_name, file_path, file_hash, import_date)
     VALUES (?, ?, ?, ?, ?, ?)",
    params = list(project_id, file_type, file_name, file_path, file_hash, current_time)
  )

  # Get the ID of the newly registered file
  file_id <- DBI::dbGetQuery(
    con,
    "SELECT file_id FROM input_files
     WHERE project_id = ? AND file_name = ? AND file_path = ?",
    params = list(project_id, file_name, file_path)
  )$file_id[1]

  if (verbose) message("Registered file '", file_name, "' with ID ", file_id)

  return(file_id)
}

#' List all input files for a project
#'
#' @param con A database connection object.
#' @param project_id The ID of the project.
#'
#' @return A data frame containing information about all input files for the project.
#'
#' @importFrom DBI dbGetQuery
#' @export
list_input_files <- function(con, project_id) {
  # Check if project exists
  project <- DBI::dbGetQuery(
    con,
    "SELECT project_id FROM projects WHERE project_id = ?",
    params = list(project_id)
  )

  if (nrow(project) == 0) {
    stop("Project with ID ", project_id, " not found.")
  }

  # Get files
  DBI::dbGetQuery(
    con,
    "SELECT * FROM input_files WHERE project_id = ? ORDER BY file_id",
    params = list(project_id)
  )
}

#' Get information about a specific file
#'
#' @param con A database connection object.
#' @param file_id The ID of the file.
#'
#' @return A data frame containing information about the file.
#'
#' @importFrom DBI dbGetQuery
#' @export
get_input_file <- function(con, file_id) {
  file_info <- DBI::dbGetQuery(
    con,
    "SELECT * FROM input_files WHERE file_id = ?",
    params = list(file_id)
  )

  if (nrow(file_info) == 0) {
    stop("File with ID ", file_id, " not found.")
  }

  return(file_info)
}
