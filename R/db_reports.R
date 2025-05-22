#' Create Dynamic Analysis Report
#'
#' Creates a dynamic Rmd/Qmd report that automatically updates as the database grows.
#' This report provides a living document of your analysis progress, with sections
#' for project overview, current status, and timestamped updates.
#'
#' @param con A database connection object created with connect_funseq_db()
#' @param project_id Integer. The project ID for which to create the report
#' @param report_path Character. Path where the report should be created. 
#'   If NULL (default), uses sanitized project name + "_analysis_report.rmd/qmd"
#' @param format Character. Report format: "Rmd" or "Qmd". Default is "Rmd"
#'   - "Rmd": R Markdown format, rendered with rmarkdown package
#'   - "Qmd": Quarto Markdown format, rendered with quarto package
#' @param template Character. Template complexity: "standard", "detailed", or "minimal". 
#'   Default is "standard"
#'   - "minimal": Simple table counts and basic progress
#'   - "standard": Project overview, input files, BLAST summaries, and updates
#'   - "detailed": All standard content plus plots, detailed statistics, and advanced analysis
#' @param open_report Logical. If TRUE, attempt to open the report in RStudio. Default is TRUE
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return Character. The path to the created report file
#'
#' @details
#' The report includes several key sections:
#' \itemize{
#'   \item Project Overview: Basic project information and metadata
#'   \item Current Database Status: Live table counts and statistics
#'   \item Project Summary: Input files, BLAST searches, and annotations
#'   \item Analysis Updates: Timestamped log of analysis progress (auto-updated)
#' }
#' 
#' The report is designed to be kept open in RStudio during analysis, providing
#' real-time insight into analysis progress and results.
#'
#' @examples
#' \dontrun{
#' # Basic usage - create standard Rmd report
#' con <- connect_funseq_db("my_analysis.db")
#' project_id <- 1
#' 
#' report_path <- create_analysis_report(con, project_id)
#' # Creates: "MyProject_analysis_report.Rmd"
#' 
#' # Create Quarto report with detailed template
#' create_analysis_report(
#'   con, project_id, 
#'   format = "Qmd", 
#'   template = "detailed"
#' )
#' # Creates: "MyProject_analysis_report.qmd" with advanced content
#' 
#' # Specify custom path and minimal template
#' create_analysis_report(
#'   con, project_id,
#'   report_path = "reports/snapper_analysis.Rmd",
#'   template = "minimal",
#'   open_report = FALSE  # Don't auto-open
#' )
#' 
#' # Create report for multiple projects
#' for (pid in c(1, 2, 3)) {
#'   create_analysis_report(con, pid, template = "standard")
#' }
#' }
#'
#' @seealso 
#' \code{\link{refresh_analysis_report}} for updating existing reports,
#' \code{\link{update_analysis_report}} for programmatic updates,
#' \code{\link{get_project_summary}} for project information
#'
#' @export
create_analysis_report <- function(con, project_id, report_path = NULL, format = "Rmd", 
                                   template = "standard", open_report = TRUE, verbose = TRUE) {
  
  # Get project info
  project_info <- get_project(con, project_id)
  
  # Generate report path if not provided
  if (is.null(report_path)) {
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", project_info$project_name)
    report_path <- paste0(safe_name, "_analysis_report.", tolower(format))
  }
  
  if (verbose) message("Creating dynamic analysis report: ", report_path)
  
  # Store report info in database for future updates
  store_report_info(con, project_id, report_path, format, template)
  
  # Create the report template
  template_content <- generate_report_template(project_info, format, template)
  
  # Write the template
  writeLines(template_content, report_path)
  
  # Initialize with current database summary
  update_analysis_report(con, project_id, section = "initialization", 
                         message = "Analysis report created", verbose = FALSE)
  
  if (verbose) message("Report created successfully: ", report_path)
  
  # Open in RStudio if requested and available
  if (open_report && interactive()) {
    tryCatch({
      if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
        rstudioapi::navigateToFile(report_path)
        if (verbose) message("Opened report in RStudio")
      } else {
        if (verbose) message("RStudio API not available - please open report manually")
      }
    }, error = function(e) {
      if (verbose) message("Could not open report automatically: ", e$message)
    })
  }
  
  return(report_path)
}

#' Update Analysis Report
#'
#' Adds a new section or updates existing sections in the dynamic report
#'
#' @param con A database connection object
#' @param project_id The project ID
#' @param section The section name (e.g., "vcf_import", "blast_search", "annotation")
#' @param message Custom message to include
#' @param auto_summary Include automatic database summary. Default is TRUE
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @export
update_analysis_report <- function(con, project_id, section = "update", message = NULL, 
                                   auto_summary = TRUE, verbose = TRUE) {
  
  # Get report info from database
  report_info <- get_report_info(con, project_id)
  
  if (is.null(report_info) || !file.exists(report_info$report_path)) {
    if (verbose) message("No report found for project ", project_id, ". Create one with create_analysis_report()")
    return(invisible(NULL))
  }
  
  if (verbose) message("Updating analysis report: ", section)
  
  # Read current report
  report_lines <- readLines(report_info$report_path, warn = FALSE)
  
  # Find the updates section or create it
  update_marker <- "<!-- DYNAMIC_UPDATES_START -->"
  update_end <- "<!-- DYNAMIC_UPDATES_END -->"
  
  update_start_idx <- which(grepl(update_marker, report_lines, fixed = TRUE))
  update_end_idx <- which(grepl(update_end, report_lines, fixed = TRUE))
  
  # Generate new update content
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  new_content <- generate_update_content(con, project_id, section, message, 
                                         timestamp, auto_summary, report_info$format)
  
  if (length(update_start_idx) > 0 && length(update_end_idx) > 0) {
    # Insert into existing updates section
    before_updates <- report_lines[1:update_start_idx]
    after_updates <- report_lines[update_end_idx:length(report_lines)]
    
    # Get existing updates
    existing_updates <- report_lines[(update_start_idx + 1):(update_end_idx - 1)]
    
    # Combine old and new updates
    updated_lines <- c(before_updates, existing_updates, new_content, after_updates)
  } else {
    # Add updates section at the end
    updated_lines <- c(
      report_lines,
      "",
      update_marker,
      new_content,
      update_end
    )
  }
  
  # Write updated report
  writeLines(updated_lines, report_info$report_path)
  
  if (verbose) message("Report updated successfully")
  return(invisible(report_info$report_path))
}

#' Generate Report Template
#'
#' @param project_info Project information
#' @param format Report format ("Rmd" or "Qmd")
#' @param template Template type
#' @return Character vector of template content
generate_report_template <- function(project_info, format, template) {
  
  # YAML header
  if (format == "Qmd") {
    yaml_header <- c(
      "---",
      paste0("title: \"Analysis Report: ", project_info$project_name, "\""),
      paste0("subtitle: \"Project ID: ", project_info$project_id, "\""),
      paste0("date: \"", format(Sys.time(), "%Y-%m-%d"), "\""),
      "format:",
      "  html:",
      "    toc: true",
      "    toc-depth: 3",
      "    code-fold: true",
      "    theme: cosmo",
      "execute:",
      "  echo: false",
      "  warning: false",
      "  message: false",
      "---"
    )
  } else {
    yaml_header <- c(
      "---",
      paste0("title: \"Analysis Report: ", project_info$project_name, "\""),
      paste0("subtitle: \"Project ID: ", project_info$project_id, "\""),
      paste0("date: \"", format(Sys.time(), "%Y-%m-%d"), "\""),
      "output:",
      "  html_document:",
      "    toc: true",
      "    toc_depth: 3",
      "    code_folding: hide",
      "    theme: cosmo",
      "---"
    )
  }
  
  # R setup chunk
  setup_chunk <- c(
    "",
    "```{r setup, include=FALSE}",
    "library(funseqR)",
    "library(DBI)",
    "library(knitr)",
    "library(ggplot2)",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
    "",
    "# Connect to database",
    paste0("con <- connect_funseq_db('", DBI::dbGetInfo(con)$dbname, "')"),
    paste0("project_id <- ", project_info$project_id),
    "```",
    ""
  )
  
  # Content based on template
  if (template == "detailed") {
    content <- generate_detailed_template(project_info)
  } else if (template == "minimal") {
    content <- generate_minimal_template(project_info)
  } else {
    content <- generate_standard_template(project_info)
  }
  
  # Updates section
  updates_section <- c(
    "",
    "# Analysis Updates",
    "",
    "This section is automatically updated as the analysis progresses.",
    "",
    "<!-- DYNAMIC_UPDATES_START -->",
    "<!-- DYNAMIC_UPDATES_END -->",
    ""
  )
  
  return(c(yaml_header, setup_chunk, content, updates_section))
}

#' Generate Standard Template Content
#' @param project_info Project information
#' @return Character vector of template content
generate_standard_template <- function(project_info) {
  c(
    "# Project Overview",
    "",
    paste0("**Project Name:** ", project_info$project_name),
    paste0("**Description:** ", project_info$description %||% "No description provided"),
    paste0("**Created:** ", project_info$creation_date),
    paste0("**Last Modified:** ", project_info$last_modified),
    "",
    "# Current Database Status",
    "",
    "```{r database-summary}",
    "summary <- get_database_summary(con, verbose = FALSE)",
    "kable(data.frame(",
    "  Table = names(summary$table_counts),",
    "  Records = unlist(summary$table_counts)",
    "), caption = 'Database Table Counts')",
    "```",
    "",
    "# Project Summary",
    "",
    "```{r project-summary}",
    "proj_summary <- get_project_summary(con, project_id, verbose = FALSE)",
    "```",
    "",
    "## Input Files",
    "",
    "```{r input-files}",
    "if (nrow(proj_summary$input_files) > 0) {",
    "  kable(proj_summary$input_files[, c('file_name', 'file_type', 'import_date')],",
    "        caption = 'Input Files')",
    "} else {",
    "  cat('No input files yet.')",
    "}",
    "```",
    "",
    "## BLAST Analyses",
    "",
    "```{r blast-summary}",
    "if (nrow(proj_summary$blast_summary) > 0) {",
    "  blast_display <- proj_summary$blast_summary[, c('blast_param_id', 'blast_type', 'db_name', 'execution_date', 'result_count')]",
    "  kable(blast_display, caption = 'BLAST Searches')",
    "} else {",
    "  cat('No BLAST searches yet.')",
    "}",
    "```",
    ""
  )
}

#' Generate Update Content
#' @param con Database connection
#' @param project_id Project ID
#' @param section Section name
#' @param message Custom message
#' @param timestamp Timestamp
#' @param auto_summary Include auto summary
#' @param format Report format
#' @return Character vector of update content
generate_update_content <- function(con, project_id, section, message, timestamp, auto_summary, format) {
  
  content <- c(
    "",
    paste0("## ", format(as.POSIXct(timestamp), "%Y-%m-%d %H:%M:%S"), " - ", stringr::str_to_title(gsub("_", " ", section)))
  )
  
  if (!is.null(message)) {
    content <- c(content, "", message, "")
  }
  
  if (auto_summary) {
    content <- c(content,
                 "",
                 "```{r}",
                 paste0("# Auto-generated summary for ", section),
                 "summary <- get_project_summary(con, project_id, verbose = FALSE)",
                 "cat('VCF Entries:', sum(summary$vcf_summary$variant_count, na.rm = TRUE), '\\n')",
                 "cat('BLAST Results:', sum(summary$blast_summary$result_count, na.rm = TRUE), '\\n')",
                 "cat('Annotations:', sum(summary$blast_summary$annotation_count, na.rm = TRUE), '\\n')",
                 "```",
                 ""
    )
  }
  
  return(content)
}

#' Store Report Information in Database
#' @param con Database connection
#' @param project_id Project ID
#' @param report_path Report path
#' @param format Report format
#' @param template Template type
store_report_info <- function(con, project_id, report_path, format, template) {
  
  # Ensure reports table exists
  ensure_reports_table(con)
  
  # Check if report already exists
  existing <- DBI::dbGetQuery(con,
                              "SELECT report_id FROM analysis_reports WHERE project_id = ?",
                              params = list(project_id))
  
  current_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  if (nrow(existing) > 0) {
    # Update existing
    DBI::dbExecute(con,
                   "UPDATE analysis_reports SET report_path = ?, format = ?, template = ?, last_updated = ? WHERE project_id = ?",
                   params = list(report_path, format, template, current_time, project_id))
  } else {
    # Insert new
    DBI::dbExecute(con,
                   "INSERT INTO analysis_reports (project_id, report_path, format, template, created_date, last_updated) VALUES (?, ?, ?, ?, ?, ?)",
                   params = list(project_id, report_path, format, template, current_time, current_time))
  }
}

#' Get Report Information from Database
#' @param con Database connection
#' @param project_id Project ID
#' @return List with report information or NULL
get_report_info <- function(con, project_id) {
  
  # Check if table exists
  tables <- DBI::dbListTables(con)
  if (!"analysis_reports" %in% tables) {
    return(NULL)
  }
  
  report_info <- DBI::dbGetQuery(con,
                                 "SELECT * FROM analysis_reports WHERE project_id = ?",
                                 params = list(project_id))
  
  if (nrow(report_info) == 0) {
    return(NULL)
  }
  
  return(as.list(report_info[1, ]))
}

#' Ensure Reports Table Exists
#' @param con Database connection
ensure_reports_table <- function(con) {
  tables <- DBI::dbListTables(con)
  if ("analysis_reports" %in% tables) {
    return(invisible(NULL))
  }
  
  DBI::dbExecute(con, "
    CREATE TABLE analysis_reports (
      report_id INTEGER PRIMARY KEY,
      project_id INTEGER NOT NULL,
      report_path TEXT NOT NULL,
      format TEXT NOT NULL,
      template TEXT NOT NULL,
      created_date TEXT NOT NULL,
      last_updated TEXT NOT NULL,
      FOREIGN KEY (project_id) REFERENCES projects (project_id)
    )
  ")
  
  DBI::dbExecute(con, "CREATE INDEX idx_reports_project ON analysis_reports (project_id)")
}

# Minimal and detailed template generators (simplified versions)
generate_minimal_template <- function(project_info) {
  c(
    "# Analysis Progress",
    "",
    paste0("Project: **", project_info$project_name, "**"),
    "",
    "```{r}",
    "get_table_counts(con, verbose = FALSE)",
    "```"
  )
}

generate_detailed_template <- function(project_info) {
  # This would be much longer with detailed plots, statistics, etc.
  c(generate_standard_template(project_info),
    "",
    "# Detailed Analysis",
    "",
    "```{r plots}",
    "# Add plots here when data is available",
    "```")
}

#' Refresh Analysis Report
#'
#' Updates an existing analysis report with current database status and optionally
#' renders it to HTML. This function is useful for getting an up-to-date view of
#' your analysis progress, especially when run interactively.
#'
#' @param con A database connection object created with connect_funseq_db()
#' @param project_id Integer. The project ID whose report should be refreshed
#' @param render_html Logical. If TRUE, render the report to HTML format after updating.
#'   Default is FALSE. Requires rmarkdown (for .Rmd) or quarto (for .Qmd) packages
#' @param open_html Logical. If TRUE and render_html=TRUE, automatically open the 
#'   rendered HTML file in the default web browser. Default is TRUE
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return Invisibly returns the path to the report file (character) or NULL if no report exists
#'
#' @details
#' This function performs several actions:
#' \itemize{
#'   \item Adds a new "Manual Refresh" section to the report with current timestamp
#'   \item Updates all dynamic content (table counts, project summaries, etc.)
#'   \item Optionally renders the updated report to HTML for easy viewing
#'   \item Can automatically open the HTML version in your web browser
#' }
#' 
#' The function will fail gracefully if no report exists for the specified project,
#' suggesting the user create one with \code{create_analysis_report()}.
#'
#' @section Rendering Requirements:
#' To use \code{render_html = TRUE}, you need:
#' \itemize{
#'   \item For .Rmd files: \code{rmarkdown} package installed
#'   \item For .Qmd files: \code{quarto} package and Quarto CLI installed
#'   \item All packages used in the report (funseqR, DBI, ggplot2, knitr, etc.)
#' }
#'
#' @examples
#' \dontrun{
#' # Connect to database
#' con <- connect_funseq_db("my_analysis.db")
#' project_id <- 1
#' 
#' # Simple refresh - just update the report content
#' refresh_analysis_report(con, project_id)
#' # Adds new section: "2025-05-22 16:23:45 - Manual Refresh"
#' 
#' # Refresh and render to HTML for presentation
#' refresh_analysis_report(con, project_id, render_html = TRUE)
#' # Updates report AND creates HTML file, opens in browser
#' 
#' # Refresh and render but don't auto-open browser
#' refresh_analysis_report(
#'   con, project_id, 
#'   render_html = TRUE, 
#'   open_html = FALSE
#' )
#' 
#' # Batch refresh multiple projects
#' project_ids <- c(1, 2, 3)
#' lapply(project_ids, function(pid) {
#'   refresh_analysis_report(con, pid, render_html = TRUE)
#' })
#' 
#' # Silent refresh (no console output)
#' refresh_analysis_report(con, project_id, verbose = FALSE)
#' 
#' # Use in a monitoring loop
#' while (analysis_running) {
#'   Sys.sleep(300)  # Wait 5 minutes
#'   refresh_analysis_report(con, project_id, render_html = TRUE)
#' }
#' }
#'
#' @section Common Use Cases:
#' \itemize{
#'   \item \strong{Progress Monitoring:} Run periodically during long analyses
#'   \item \strong{Presentation Prep:} Render to HTML before meetings/presentations
#'   \item \strong{Checkpoint Documentation:} Create timestamped snapshots of progress
#'   \item \strong{Troubleshooting:} Get current status when debugging issues
#' }
#'
#' @seealso 
#' \code{\link{create_analysis_report}} for creating new reports,
#' \code{\link{update_analysis_report}} for programmatic updates,
#' \code{\link{get_database_summary}} for raw database statistics
#'
#' @export
refresh_analysis_report <- function(con, project_id, render_html = FALSE, open_html = TRUE) {
  
  # Update with current status
  update_analysis_report(
    con, project_id,
    section = "manual_refresh",
    message = "Report manually refreshed with current database status",
    auto_summary = TRUE,
    verbose = TRUE
  )
  
  # Render to HTML if requested
  if (render_html) {
    report_info <- get_report_info(con, project_id)
    
    if (!is.null(report_info) && file.exists(report_info$report_path)) {
      tryCatch({
        if (report_info$format == "Qmd" && requireNamespace("quarto", quietly = TRUE)) {
          quarto::quarto_render(report_info$report_path)
          html_file <- gsub("\\.qmd$", ".html", report_info$report_path)
        } else if (requireNamespace("rmarkdown", quietly = TRUE)) {
          rmarkdown::render(report_info$report_path, quiet = TRUE)
          html_file <- gsub("\\.(Rmd|rmd)$", ".html", report_info$report_path)
        }
        
        if (open_html && exists("html_file") && file.exists(html_file)) {
          browseURL(html_file)
          message("Opened rendered report: ", html_file)
        }
        
      }, error = function(e) {
        message("Could not render report: ", e$message)
      })
    }
  }
}