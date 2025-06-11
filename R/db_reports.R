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
  template_content <- generate_report_template(con, project_info, format, template)
  
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
    # Use the first valid pair of markers
    start_idx <- update_start_idx[1]
    end_idx <- update_end_idx[length(update_end_idx)]  # Use last end marker
    
    # Insert into existing updates section
    before_updates <- report_lines[1:start_idx]
    after_updates <- report_lines[end_idx:length(report_lines)]
    
    # Get existing updates (between first start and last end)
    if (end_idx > start_idx + 1) {
      existing_updates <- report_lines[(start_idx + 1):(end_idx - 1)]
      # Remove any duplicate markers from existing content
      existing_updates <- existing_updates[!grepl("<!-- DYNAMIC_UPDATES_(START|END) -->", existing_updates)]
    } else {
      existing_updates <- character(0)
    }
    
    # Combine old and new updates
    updated_lines <- c(before_updates, new_content, existing_updates, after_updates)
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
#' @param con Database connection
#' @param project_info Project information
#' @param format Report format ("Rmd" or "Qmd")
#' @param template Template type
#' @return Character vector of template content
generate_report_template <- function(con, project_info, format, template) {
  
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
  content <- generate_template_content(project_info, template)
  
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

#' Generate Template Content Based on Type
#' @param project_info Project information
#' @param template_type Template type ("minimal", "standard", "detailed")
#' @return Character vector of template content
generate_template_content <- function(project_info, template_type) {
  content <- character(0)
  
  # Minimal template - only basic content
  if (template_type == "minimal") {
    return(c(
      "# Analysis Progress",
      "",
      paste0("Project: **", project_info$project_name, "**"),
      "",
      "```{r}",
      "get_table_counts(con, verbose = FALSE)",
      "```"
    ))
  }
  
  # Standard and detailed templates start with comprehensive content
  if (template_type %in% c("standard", "detailed")) {
    content <- c(
      # Project Overview Section
      "# Project Overview",
      "",
      paste0("**Project Name:** ", project_info$project_name),
      paste0("**Description:** ", project_info$description %||% "No description provided"),
      paste0("**Created:** ", project_info$creation_date),
      paste0("**Last Modified:** ", project_info$last_modified),
      "",
      
      # Database Status Section
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
      
      # Workflow Status Section
      "# Analysis Workflow Status",
      "",
      "```{r workflow-status}",
      "proj_summary <- get_project_summary(con, project_id, verbose = FALSE)",
      "",
      "workflow_steps <- c('1. VCF Import', '2. Reference Genome', '3. Flanking Sequences', '4. BLAST Search', '5. Annotation')",
      "",
      "vcf_count <- tryCatch(summary$table_counts$vcf_data, error = function(e) 0)",
      "genome_count <- tryCatch(summary$table_counts$reference_sequences, error = function(e) 0)",
      "flanking_count <- tryCatch(summary$table_counts$flanking_sequences, error = function(e) 0)",
      "blast_count <- tryCatch(summary$table_counts$blast_results, error = function(e) 0)",
      "annotation_count <- tryCatch(summary$table_counts$annotations, error = function(e) 0)",
      "",
      "workflow_status_values <- c(",
      "  ifelse(vcf_count > 0, '✓ Complete', '○ Pending'),",
      "  ifelse(genome_count > 0, '✓ Complete', '○ Pending'),",
      "  ifelse(flanking_count > 0, '✓ Complete', '○ Pending'),",
      "  ifelse(blast_count > 0, '✓ Complete', '○ Pending'),",
      "  ifelse(annotation_count > 0, '✓ Complete', '○ Pending')",
      ")",
      "",
      "workflow_records_values <- c(",
      "  format(vcf_count, big.mark = ','),",
      "  format(genome_count, big.mark = ','),",
      "  format(flanking_count, big.mark = ','),",
      "  format(blast_count, big.mark = ','),",
      "  format(annotation_count, big.mark = ',')",
      ")",
      "",
      "workflow_status <- data.frame(",
      "  Step = workflow_steps,",
      "  Status = workflow_status_values,",
      "  Records = workflow_records_values,",
      "  stringsAsFactors = FALSE",
      ")",
      "",
      "kable(workflow_status, caption = 'Analysis Workflow Progress')",
      "```",
      "",
      
      # Input Data Summary Section
      "# Input Data Summary",
      "",
      "## VCF Files",
      "",
      "```{r vcf-files}",
      "if ('input_files' %in% names(proj_summary) && nrow(proj_summary$input_files) > 0) {",
      "  vcf_files <- proj_summary$input_files[tolower(proj_summary$input_files$file_type) == 'vcf', ]",
      "  if (nrow(vcf_files) > 0) {",
      "    if ('vcf_summary' %in% names(proj_summary) && 'file_id' %in% names(vcf_files) && 'file_id' %in% names(proj_summary$vcf_summary)) {",
      "      vcf_display <- merge(vcf_files, proj_summary$vcf_summary, by = 'file_id', all.x = TRUE)",
      "      vcf_display <- vcf_display[, c('file_name', 'variant_count', 'import_date')]",
      "      kable(vcf_display, caption = 'VCF Files Imported')",
      "    } else {",
      "      vcf_display <- vcf_files[, c('file_name', 'import_date')]",
      "      kable(vcf_display, caption = 'VCF Files Imported')",
      "    }",
      "  } else {",
      "    cat('No VCF files imported yet.')",
      "  }",
      "} else {",
      "  cat('No VCF files imported yet.')",
      "}",
      "```",
      "",
      "## Reference Genomes",
      "",
      "```{r reference-genomes}",
      "if ('genome_summary' %in% names(proj_summary) && nrow(proj_summary$genome_summary) > 0) {",
      "  genome_display <- proj_summary$genome_summary[, c('genome_name', 'genome_build', 'sequence_count', 'import_date')]",
      "  kable(genome_display, caption = 'Reference Genomes')",
      "} else {",
      "  cat('No reference genomes imported yet.')",
      "}",
      "```",
      "",
      
      # Analysis Results Section
      "# Analysis Results",
      "",
      "## BLAST Searches",
      "",
      "```{r blast-summary}",
      "if ('blast_summary' %in% names(proj_summary) && nrow(proj_summary$blast_summary) > 0) {",
      "  blast_display <- proj_summary$blast_summary[, c('blast_param_id', 'blast_type', 'db_name', 'execution_date', 'result_count')]",
      "  kable(blast_display, caption = 'BLAST Searches Performed')",
      "} else {",
      "  cat('No BLAST searches performed yet.')",
      "}",
      "```",
      "",
      "## Functional Annotations",
      "",
      "```{r annotation-summary}",
      "if (summary$table_counts$annotations > 0) {",
      "  annotation_stats <- DBI::dbGetQuery(con, '",
      "    SELECT ",
      "      COUNT(*) as total_annotations,",
      "      COUNT(DISTINCT uniprot_accession) as unique_proteins,",
      "      (SELECT COUNT(*) FROM go_terms) as go_terms,",
      "      (SELECT COUNT(*) FROM kegg_references) as kegg_refs",
      "    FROM annotations",
      "  ')",
      "  ",
      "  annotation_display <- data.frame(",
      "    Metric = c('Total Annotations', 'Unique Proteins', 'GO Terms', 'KEGG References'),",
      "    Count = c(",
      "      format(annotation_stats$total_annotations, big.mark = ','),",
      "      format(annotation_stats$unique_proteins, big.mark = ','),",
      "      format(annotation_stats$go_terms, big.mark = ','),",
      "      format(annotation_stats$kegg_refs, big.mark = ',')",
      "    )",
      "  )",
      "  ",
      "  kable(annotation_display, caption = 'Functional Annotation Summary')",
      "} else {",
      "  cat('No functional annotations available yet.')",
      "}",
      "```",
      ""
    )
  }
  
  # Add detailed sections for detailed template
  if (template_type == "detailed") {
    detailed_content <- c(
      "",
      "# Detailed Analysis",
      "",
      "## Chromosome Distribution",
      "",
      "```{r chromosome-distribution, fig.height=6, fig.width=10}",
      "if (summary$table_counts$vcf_data > 0) {",
      "  chrom_dist <- DBI::dbGetQuery(con, '",
      "    SELECT chromosome, COUNT(*) as variant_count",
      "    FROM vcf_data v",
      "    JOIN input_files if ON v.file_id = if.file_id",
      "    WHERE if.project_id = ?",
      "    GROUP BY chromosome",
      "    ORDER BY chromosome",
      "  ', params = list(project_id))",
      "  ",
      "  if (nrow(chrom_dist) > 0) {",
      "    p1 <- ggplot(chrom_dist, aes(x = reorder(chromosome, -variant_count), y = variant_count)) +",
      "      geom_bar(stat = 'identity', fill = 'steelblue', alpha = 0.7) +",
      "      labs(title = 'Variant Distribution by Chromosome',",
      "           x = 'Chromosome', y = 'Number of Variants') +",
      "      theme_minimal() +",
      "      theme(axis.text.x = element_text(angle = 45, hjust = 1))",
      "    print(p1)",
      "  }",
      "} else {",
      "  cat('No variant data available for chromosome distribution plot.')",
      "}",
      "```",
      "",
      "## BLAST Results Quality",
      "",
      "```{r blast-quality, fig.height=6, fig.width=10}",
      "if (summary$table_counts$blast_results > 0) {",
      "  blast_stats <- DBI::dbGetQuery(con, '",
      "    SELECT br.e_value, br.bit_score, br.percent_identity",
      "    FROM blast_results br",
      "    JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id",
      "    WHERE bp.project_id = ?",
      "  ', params = list(project_id))",
      "  ",
      "  if (nrow(blast_stats) > 0) {",
      "    p2 <- ggplot(blast_stats, aes(x = log10(e_value))) +",
      "      geom_histogram(bins = 30, fill = 'orange', alpha = 0.7) +",
      "      labs(title = 'BLAST E-value Distribution',",
      "           x = 'log10(E-value)', y = 'Frequency') +",
      "      theme_minimal()",
      "    ",
      "    p3 <- ggplot(blast_stats, aes(x = percent_identity, y = bit_score)) +",
      "      geom_point(alpha = 0.5, color = 'darkgreen') +",
      "      labs(title = 'BLAST Hit Quality: Identity vs Bit Score',",
      "           x = 'Percent Identity', y = 'Bit Score') +",
      "      theme_minimal()",
      "    ",
      "    print(p2)",
      "    print(p3)",
      "  }",
      "} else {",
      "  cat('No BLAST results available for quality plots.')",
      "}",
      "```",
      "",
      "## GO Term Analysis",
      "",
      "```{r go-analysis, fig.height=6, fig.width=10}",
      "if (summary$table_counts$go_terms > 0) {",
      "  go_dist <- DBI::dbGetQuery(con, '",
      "    SELECT go_category, COUNT(*) as term_count",
      "    FROM go_terms gt",
      "    JOIN annotations a ON gt.annotation_id = a.annotation_id",
      "    JOIN blast_results br ON a.blast_result_id = br.blast_result_id",
      "    JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id",
      "    WHERE bp.project_id = ?",
      "    GROUP BY go_category",
      "  ', params = list(project_id))",
      "  ",
      "  if (nrow(go_dist) > 0) {",
      "    go_dist$category_name <- ifelse(go_dist$go_category == 'P', 'Biological Process',",
      "                                   ifelse(go_dist$go_category == 'F', 'Molecular Function',",
      "                                         ifelse(go_dist$go_category == 'C', 'Cellular Component', 'Unknown')))",
      "    ",
      "    p4 <- ggplot(go_dist, aes(x = '', y = term_count, fill = category_name)) +",
      "      geom_bar(stat = 'identity', width = 1) +",
      "      coord_polar('y', start = 0) +",
      "      labs(title = 'GO Term Distribution by Category', fill = 'GO Category') +",
      "      theme_void()",
      "    print(p4)",
      "  }",
      "} else {",
      "  cat('No GO terms available for analysis.')",
      "}",
      "```",
      "",
      "## Data Quality Metrics",
      "",
      "```{r quality-metrics}",
      "quality_metrics <- data.frame(",
      "  Metric = character(0),",
      "  Value = character(0)",
      ")",
      "",
      "if (summary$table_counts$vcf_data > 0 && summary$table_counts$flanking_sequences > 0) {",
      "  flanking_success <- round((summary$table_counts$flanking_sequences / summary$table_counts$vcf_data) * 100, 1)",
      "  quality_metrics <- rbind(quality_metrics, data.frame(",
      "    Metric = 'Flanking Sequence Success Rate',",
      "    Value = paste0(flanking_success, '%')",
      "  ))",
      "}",
      "",
      "if (summary$table_counts$flanking_sequences > 0 && summary$table_counts$blast_results > 0) {",
      "  blast_hit_rate <- round((summary$table_counts$blast_results / summary$table_counts$flanking_sequences) * 100, 1)",
      "  quality_metrics <- rbind(quality_metrics, data.frame(",
      "    Metric = 'BLAST Hit Rate',",
      "    Value = paste0(blast_hit_rate, '%')",
      "  ))",
      "}",
      "",
      "if (summary$table_counts$blast_results > 0 && summary$table_counts$annotations > 0) {",
      "  annotation_rate <- round((summary$table_counts$annotations / summary$table_counts$blast_results) * 100, 1)",
      "  quality_metrics <- rbind(quality_metrics, data.frame(",
      "    Metric = 'Annotation Success Rate',",
      "    Value = paste0(annotation_rate, '%')",
      "  ))",
      "}",
      "",
      "if (nrow(quality_metrics) > 0) {",
      "  kable(quality_metrics, caption = 'Analysis Quality Metrics')",
      "} else {",
      "  cat('Insufficient data for quality metrics.')",
      "}",
      "```"
    )
    
    content <- c(content, detailed_content)
  }
  
  return(content)
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

#' Add GO Enrichment Summary to Dynamic Report
#' 
#' @param con Database connection
#' @param project_id Project ID
#' @param go_content_file Path to GO report content RDS file
#' @param verbose Logical. Print progress information
#' 
#' @export
add_go_enrichment_to_report <- function(con, project_id, go_content_file, verbose = TRUE) {
  
  if (verbose) message("Adding GO enrichment summary to dynamic report...")
  
  # Check if GO content file exists
  if (!file.exists(go_content_file)) {
    stop("GO content file not found: ", go_content_file)
  }
  
  # Load GO content
  go_content <- readRDS(go_content_file)
  
  # Get report info
  report_info <- get_report_info(con, project_id)
  
  if (is.null(report_info)) {
    stop("No report found for project ", project_id, ". Create a report first.")
  }
  
  # Read current report
  if (!file.exists(report_info$report_path)) {
    stop("Report file not found: ", report_info$report_path)
  }
  
  report_lines <- readLines(report_info$report_path)
  
  # Find insertion point (before "# Analysis Updates" section)
  updates_idx <- which(grepl("^# Analysis Updates", report_lines))
  
  if (length(updates_idx) == 0) {
    # Add at end if no updates section found
    insertion_point <- length(report_lines) + 1
  } else {
    insertion_point <- updates_idx[1]
  }
  
  # Create GO enrichment section
  go_section <- c(
    "",
    "# GO Enrichment Analysis",
    "",
    paste0("*Analysis Date: ", go_content$analysis_date, "*"),
    "",
    go_content$summary_text,
    ""
  )
  
  # Insert GO section
  new_report <- c(
    report_lines[1:(insertion_point-1)],
    go_section,
    report_lines[insertion_point:length(report_lines)]
  )
  
  # Write updated report
  writeLines(new_report, report_info$report_path)
  
  # Update timestamp in database
  current_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  DBI::dbExecute(con,
                 "UPDATE analysis_reports SET last_updated = ? WHERE project_id = ?",
                 list(current_time, project_id))
  
  if (verbose) {
    message("GO enrichment summary added to report: ", report_info$report_path)
    message("Report updated with ", go_content$total_enriched, " enriched GO terms")
  }
  
  return(invisible(report_info$report_path))
}