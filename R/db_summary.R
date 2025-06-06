#' Get comprehensive database summary
#'
#' This function provides a summary of the entire funseqR database contents
#'
#' @param con A database connection object
#' @param detail_level Character. Level of detail: "basic", "detailed", or "full". Default is "detailed"
#' @param include_samples Logical. Include sample data from tables. Default is FALSE
#' @param verbose Logical. Print detailed output. Default is TRUE
#'
#' @return A list containing database summary information
#' @export
get_database_summary <- function(con, detail_level = "detailed", include_samples = FALSE, verbose = TRUE) {
  if (!DBI::dbIsValid(con)) {
    stop("Invalid database connection")
  }

  # Helper function for null coalescing
  `%||%` <- function(x, y) if (is.null(x) || is.na(x)) y else x

  # Get all tables in the database
  all_tables <- DBI::dbListTables(con)

  if (verbose) {
    cat("=== funseqR Database Summary ===\n")
    cat("Tables found:", length(all_tables), "\n")
    cat("Tables:", paste(all_tables, collapse = ", "), "\n\n")
  }

  # Get basic table counts
  table_counts <- list()

  # Core tables with counts
  core_tables <- c("projects", "input_files", "vcf_data", "reference_genomes",
                   "reference_sequences", "flanking_sequences", "blast_parameters",
                   "blast_results", "annotations", "go_terms", "kegg_references",
                   "uniprot_cache", "blast_database_metadata")

  for (table in core_tables) {
    if (table %in% all_tables) {
      count_query <- paste0("SELECT COUNT(*) as count FROM ", table)
      tryCatch({
        count <- DBI::dbGetQuery(con, count_query)$count
        table_counts[[table]] <- count
        if (verbose) cat(sprintf("%-25s: %d records\n", table, count))
      }, error = function(e) {
        table_counts[[table]] <- NA
        if (verbose) cat(sprintf("%-25s: Error reading\n", table))
      })
    } else {
      table_counts[[table]] <- 0
      if (verbose) cat(sprintf("%-25s: Table not found\n", table))
    }
  }

  # Get metadata if available
  metadata <- list()
  if ("metadata" %in% all_tables) {
    tryCatch({
      metadata_query <- "SELECT * FROM metadata"
      metadata_raw <- DBI::dbGetQuery(con, metadata_query)
      metadata <- setNames(metadata_raw$value, metadata_raw$key)

      if (verbose && length(metadata) > 0) {
        cat("\n=== Database Metadata ===\n")
        for (key in names(metadata)) {
          cat(sprintf("%-20s: %s\n", key, metadata[[key]]))
        }
      }
    }, error = function(e) {
      if (verbose) cat("Could not read metadata table\n")
    })
  }

  # Get project-level summary
  project_summary <- list()
  if (table_counts$projects > 0) {
    if (verbose) cat("\n=== Project Summary ===\n")

    tryCatch({
      project_query <- "
        SELECT p.project_id, p.project_name, p.creation_date,
               COUNT(DISTINCT f.file_id) as input_files,
               COUNT(DISTINCT v.vcf_id) as vcf_entries,
               COUNT(DISTINCT bp.blast_param_id) as blast_runs,
               COUNT(DISTINCT a.annotation_id) as annotations
        FROM projects p
        LEFT JOIN input_files f ON p.project_id = f.project_id
        LEFT JOIN vcf_data v ON f.file_id = v.file_id
        LEFT JOIN blast_parameters bp ON p.project_id = bp.project_id
        LEFT JOIN blast_results br ON bp.blast_param_id = br.blast_param_id
        LEFT JOIN annotations a ON br.blast_result_id = a.blast_result_id
        GROUP BY p.project_id, p.project_name, p.creation_date
        ORDER BY p.project_id
      "

      project_summary <- DBI::dbGetQuery(con, project_query)

      if (verbose && nrow(project_summary) > 0) {
        for (i in 1:nrow(project_summary)) {
          proj <- project_summary[i, ]
          cat(sprintf("Project %d: %s\n", proj$project_id, proj$project_name))
          cat(sprintf("  Created: %s\n", proj$creation_date))
          cat(sprintf("  Input files: %d, VCF entries: %d\n", proj$input_files, proj$vcf_entries))
          cat(sprintf("  BLAST runs: %d, Annotations: %d\n", proj$blast_runs, proj$annotations))
        }
      }

    }, error = function(e) {
      if (verbose) cat("Could not generate project summary\n")
    })
  }

  # Detailed analysis if requested
  detailed_info <- list()
  if (detail_level %in% c("detailed", "full")) {

    # GO term analysis
    if (table_counts$go_terms > 0) {
      tryCatch({
        go_summary_query <- "
          SELECT go_category, COUNT(*) as count
          FROM go_terms
          GROUP BY go_category
          ORDER BY count DESC
        "
        go_summary <- DBI::dbGetQuery(con, go_summary_query)
        detailed_info$go_categories <- go_summary

        if (verbose) {
          cat("\n=== GO Terms by Category ===\n")
          for (i in 1:nrow(go_summary)) {
            cat_name <- switch(go_summary$go_category[i],
                               "P" = "Biological Process",
                               "F" = "Molecular Function",
                               "C" = "Cellular Component",
                               go_summary$go_category[i])
            cat(sprintf("%-20s: %d\n", cat_name, go_summary$count[i]))
          }
        }
      }, error = function(e) {
        if (verbose) cat("Could not analyze GO terms\n")
      })
    }

    # BLAST database analysis
    if (table_counts$blast_parameters > 0) {
      tryCatch({
        # Basic BLAST parameters summary
        blast_db_query <- "
          SELECT bp.db_name, bp.blast_type, COUNT(*) as runs,
                 AVG(bp.e_value) as avg_evalue, MIN(bp.execution_date) as first_run,
                 MAX(bp.execution_date) as last_run
          FROM blast_parameters bp
          GROUP BY bp.db_name, bp.blast_type
          ORDER BY runs DESC
        "
        blast_summary <- DBI::dbGetQuery(con, blast_db_query)
        detailed_info$blast_databases <- blast_summary

        # Enhanced: Get database metadata if available
        if ("blast_database_metadata" %in% all_tables && table_counts$blast_database_metadata > 0) {
          db_metadata_query <- "
            SELECT DISTINCT bm.db_name, bm.db_title, bm.num_sequences,
                   bm.total_length, bm.db_date, bm.db_version,
                   bm.extraction_date, bp.blast_type,
                   COUNT(DISTINCT bp.blast_param_id) as times_used
            FROM blast_database_metadata bm
            JOIN blast_parameters bp ON bm.blast_param_id = bp.blast_param_id
            GROUP BY bm.db_name, bm.db_title, bm.num_sequences,
                     bm.total_length, bm.db_date, bm.db_version, bp.blast_type
            ORDER BY times_used DESC, bm.extraction_date DESC
          "

          db_metadata_summary <- DBI::dbGetQuery(con, db_metadata_query)
          detailed_info$blast_database_metadata <- db_metadata_summary

          if (verbose && nrow(db_metadata_summary) > 0) {
            cat("\n=== BLAST Database Details (with Metadata) ===\n")
            for (i in 1:nrow(db_metadata_summary)) {
              bm <- db_metadata_summary[i, ]
              cat(sprintf("Database: %s (%s)\n",
                          bm$db_name, bm$blast_type))
              cat(sprintf("  Title: %s\n",
                          bm$db_title %||% "Unknown"))
              cat(sprintf("  Sequences: %s\n",
                          format(bm$num_sequences %||% 0, big.mark = ",")))
              cat(sprintf("  Total length: %s bp\n",
                          format(bm$total_length %||% 0, big.mark = ",")))
              cat(sprintf("  Database date: %s\n",
                          bm$db_date %||% "Unknown"))
              if (!is.na(bm$db_version) && !is.null(bm$db_version)) {
                cat(sprintf("  Version: %s\n", bm$db_version))
              }
              cat(sprintf("  Times used: %d\n", bm$times_used))
              cat(sprintf("  Metadata extracted: %s\n", bm$extraction_date))
              cat("\n")
            }
          }
        } else {
          # Fallback to basic summary if no metadata available
          if (verbose && nrow(blast_summary) > 0) {
            cat("\n=== BLAST Database Usage ===\n")
            for (i in 1:nrow(blast_summary)) {
              bs <- blast_summary[i, ]
              cat(sprintf("Database: %s (%s)\n", bs$db_name, bs$blast_type))
              cat(sprintf("  Runs: %d, Avg E-value: %.2e\n", bs$runs, bs$avg_evalue))
              cat(sprintf("  First run: %s, Last run: %s\n", bs$first_run, bs$last_run))
            }
          }
        }
      }, error = function(e) {
        if (verbose) cat("Could not analyze BLAST parameters\n")
      })
    }

    # Chromosome distribution (FIXED for SQLite)
    if (table_counts$vcf_data > 0) {
      tryCatch({
        # SQLite-compatible query (no REGEXP, simpler sorting)
        chrom_query <- "
          SELECT chromosome, COUNT(*) as variant_count
          FROM vcf_data
          GROUP BY chromosome
          ORDER BY
            CASE
              WHEN chromosome GLOB '[0-9]*' THEN CAST(chromosome AS INTEGER)
              WHEN chromosome LIKE 'chr%' AND SUBSTR(chromosome, 4) GLOB '[0-9]*'
                THEN CAST(SUBSTR(chromosome, 4) AS INTEGER)
              ELSE 999
            END,
            chromosome
        "
        chrom_summary <- DBI::dbGetQuery(con, chrom_query)
        detailed_info$chromosome_distribution <- chrom_summary

        if (verbose && nrow(chrom_summary) > 0) {
          cat("\n=== Variants by Chromosome ===\n")
          total_variants <- sum(chrom_summary$variant_count)
          for (i in 1:min(nrow(chrom_summary), 10)) {  # Show top 10
            cs <- chrom_summary[i, ]
            pct <- round(100 * cs$variant_count / total_variants, 1)
            cat(sprintf("%-10s: %d variants (%s%%)\n", cs$chromosome, cs$variant_count, pct))
          }
          if (nrow(chrom_summary) > 10) {
            cat(sprintf("... and %d more chromosomes\n", nrow(chrom_summary) - 10))
          }
        }
      }, error = function(e) {
        if (verbose) cat("Could not analyze chromosome distribution: ", e$message, "\n")
      })
    }
  }

  # Full analysis with sample data
  sample_data <- list()
  if (detail_level == "full" && include_samples) {

    if (verbose) cat("\n=== Sample Data ===\n")

    # Sample from each major table
    sample_tables <- c("projects", "vcf_data", "annotations", "go_terms")

    for (table in sample_tables) {
      if (table %in% all_tables && table_counts[[table]] > 0) {
        tryCatch({
          sample_query <- paste0("SELECT * FROM ", table, " LIMIT 3")
          sample_data[[table]] <- DBI::dbGetQuery(con, sample_query)

          if (verbose) {
            cat(sprintf("\n--- Sample from %s ---\n", table))
            print(head(sample_data[[table]], 3))
          }
        }, error = function(e) {
          if (verbose) cat(sprintf("Could not sample from %s\n", table))
        })
      }
    }
  }

  # Return comprehensive summary
  list(
    database_info = list(
      tables = all_tables,
      metadata = metadata
    ),
    table_counts = table_counts,
    project_summary = project_summary,
    detailed_info = detailed_info,
    sample_data = sample_data,
    summary_date = Sys.time()
  )
}

#' Get project-specific summary
#'
#' @param con A database connection object
#' @param project_id Integer. Specific project ID to analyze
#' @param verbose Logical. Print detailed output. Default is TRUE
#'
#' @return A list containing project-specific information
#' @export
get_project_summary <- function(con, project_id, verbose = TRUE) {
  if (!DBI::dbIsValid(con)) {
    stop("Invalid database connection")
  }

  # Helper function for null coalescing
  `%||%` <- function(x, y) if (is.null(x) || is.na(x)) y else x

  # Check if project exists
  project_check <- DBI::dbGetQuery(con,
                                   "SELECT * FROM projects WHERE project_id = ?",
                                   params = list(project_id))

  if (nrow(project_check) == 0) {
    stop("Project with ID ", project_id, " not found")
  }

  project_info <- project_check[1, ]

  if (verbose) {
    cat("=== Project Summary ===\n")
    cat("Project ID:", project_info$project_id, "\n")
    cat("Name:", project_info$project_name, "\n")
    cat("Description:", ifelse(is.na(project_info$description), "None", project_info$description), "\n")
    cat("Created:", project_info$creation_date, "\n")
    cat("Last modified:", project_info$last_modified, "\n\n")
  }

  # Get input files
  files_query <- "SELECT * FROM input_files WHERE project_id = ? ORDER BY file_id"
  input_files <- DBI::dbGetQuery(con, files_query, params = list(project_id))

  if (verbose && nrow(input_files) > 0) {
    cat("=== Input Files ===\n")
    for (i in 1:nrow(input_files)) {
      file <- input_files[i, ]
      cat(sprintf("File %d: %s (%s)\n", file$file_id, file$file_name, file$file_type))
      cat(sprintf("  Path: %s\n", file$file_path))
      cat(sprintf("  Import date: %s\n", file$import_date))
    }
    cat("\n")
  }

  # Get VCF data summary
  vcf_summary <- NULL
  if (nrow(input_files) > 0) {
    file_ids <- paste(input_files$file_id, collapse = ",")
    vcf_query <- paste0("
      SELECT f.file_name, COUNT(v.vcf_id) as variant_count,
             COUNT(DISTINCT v.chromosome) as chromosome_count
      FROM input_files f
      LEFT JOIN vcf_data v ON f.file_id = v.file_id
      WHERE f.project_id = ", project_id, "
      GROUP BY f.file_id, f.file_name
    ")
    vcf_summary <- DBI::dbGetQuery(con, vcf_query)

    if (verbose && nrow(vcf_summary) > 0) {
      cat("=== VCF Data Summary ===\n")
      total_variants <- sum(vcf_summary$variant_count)
      for (i in 1:nrow(vcf_summary)) {
        vs <- vcf_summary[i, ]
        cat(sprintf("%s: %d variants across %d chromosomes\n",
                    vs$file_name, vs$variant_count, vs$chromosome_count))
      }
      cat(sprintf("Total variants: %d\n\n", total_variants))
    }
  }

  # Get BLAST summary (ENHANCED)
  blast_query <- "
    SELECT bp.blast_param_id, bp.blast_type, bp.db_name, bp.execution_date,
           COUNT(br.blast_result_id) as result_count,
           COUNT(DISTINCT a.annotation_id) as annotation_count,
           bm.db_title, bm.num_sequences, bm.total_length, bm.db_date, bm.db_version
    FROM blast_parameters bp
    LEFT JOIN blast_results br ON bp.blast_param_id = br.blast_param_id
    LEFT JOIN annotations a ON br.blast_result_id = a.blast_result_id
    LEFT JOIN blast_database_metadata bm ON bp.blast_param_id = bm.blast_param_id
    WHERE bp.project_id = ?
    GROUP BY bp.blast_param_id
    ORDER BY bp.execution_date DESC
  "
  blast_summary <- DBI::dbGetQuery(con, blast_query, params = list(project_id))

  if (verbose && nrow(blast_summary) > 0) {
    cat("=== BLAST Analysis Summary ===\n")
    for (i in 1:nrow(blast_summary)) {
      bs <- blast_summary[i, ]
      cat(sprintf("BLAST Run %d (%s):\n", bs$blast_param_id, bs$execution_date))
      cat(sprintf("  Database: %s, Type: %s\n", bs$db_name, bs$blast_type))

      # Show database metadata if available
      if (!is.na(bs$db_title) && !is.null(bs$db_title)) {
        cat(sprintf("  DB Title: %s\n", bs$db_title))
        cat(sprintf("  DB Sequences: %s\n",
                    format(bs$num_sequences %||% 0, big.mark = ",")))
        cat(sprintf("  DB Length: %s bp\n",
                    format(bs$total_length %||% 0, big.mark = ",")))
        cat(sprintf("  DB Date: %s\n", bs$db_date %||% "Unknown"))
        if (!is.na(bs$db_version) && !is.null(bs$db_version)) {
          cat(sprintf("  DB Version: %s\n", bs$db_version))
        }
      }

      cat(sprintf("  Results: %d, Annotations: %d\n", bs$result_count, bs$annotation_count))
    }
    cat("\n")
  }

  # Get annotation summary
  annotation_summary <- data.frame()
  if (nrow(blast_summary) > 0) {
    annotation_query <- paste0("
      SELECT COUNT(DISTINCT a.annotation_id) as total_annotations,
             COUNT(DISTINCT g.go_term_id) as go_terms,
             COUNT(DISTINCT k.kegg_ref_id) as kegg_refs,
             COUNT(DISTINCT a.uniprot_accession) as unique_proteins
      FROM annotations a
      JOIN blast_results br ON a.blast_result_id = br.blast_result_id
      JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
      LEFT JOIN go_terms g ON a.annotation_id = g.annotation_id
      LEFT JOIN kegg_references k ON a.annotation_id = k.annotation_id
      WHERE bp.project_id = ", project_id
    )

    annotation_summary <- DBI::dbGetQuery(con, annotation_query)

    if (verbose && nrow(annotation_summary) > 0) {
      cat("=== Annotation Summary ===\n")
      as <- annotation_summary[1, ]
      cat(sprintf("Total annotations: %d\n", as$total_annotations))
      cat(sprintf("Unique proteins: %d\n", as$unique_proteins))
      cat(sprintf("GO terms: %d\n", as$go_terms))
      cat(sprintf("KEGG references: %d\n", as$kegg_refs))
    }
  }

  # Get reference genome summary
  genome_summary <- data.frame()
  tryCatch({
    genome_query <- "
      SELECT rg.genome_id, rg.genome_name, rg.genome_build,
             COUNT(DISTINCT rs.sequence_id) as sequence_count,
             if.import_date
      FROM reference_genomes rg
      JOIN input_files if ON rg.file_id = if.file_id
      LEFT JOIN reference_sequences rs ON rg.genome_id = rs.genome_id
      WHERE if.project_id = ?
      GROUP BY rg.genome_id, rg.genome_name, rg.genome_build, if.import_date
      ORDER BY if.import_date DESC
    "
    
    genome_summary <- DBI::dbGetQuery(con, genome_query, params = list(project_id))
    
    if (verbose && nrow(genome_summary) > 0) {
      cat("\n=== Reference Genomes ===\n")
      for (i in 1:nrow(genome_summary)) {
        gs <- genome_summary[i, ]
        cat(sprintf("Genome %d: %s\n", gs$genome_id, gs$genome_name))
        if (!is.na(gs$genome_build) && !is.null(gs$genome_build)) {
          cat(sprintf("  Build: %s\n", gs$genome_build))
        }
        cat(sprintf("  Sequences: %d\n", gs$sequence_count))
        cat(sprintf("  Import date: %s\n", gs$import_date))
      }
    }
  }, error = function(e) {
    if (verbose) cat("Could not retrieve reference genome information\n")
  })

  # Return comprehensive project summary
  list(
    project_info = project_info,
    input_files = input_files,
    vcf_summary = vcf_summary,
    blast_summary = blast_summary,
    annotation_summary = annotation_summary,
    genome_summary = genome_summary,
    analysis_date = Sys.time()
  )
}

#' Get detailed information about a specific table
#'
#' @param con A database connection object
#' @param table_name Character. Name of the table to analyze
#' @param sample_size Integer. Number of sample rows to return. Default is 5
#' @param include_schema Logical. Include table schema information. Default is TRUE
#' @param verbose Logical. Print detailed output. Default is TRUE
#'
#' @return A list containing table information
#' @export
get_table_info <- function(con, table_name, sample_size = 5, include_schema = TRUE, verbose = TRUE) {

  if (!DBI::dbIsValid(con)) {
    stop("Invalid database connection")
  }

  # Check if table exists
  all_tables <- DBI::dbListTables(con)
  if (!table_name %in% all_tables) {
    stop("Table '", table_name, "' not found in database")
  }

  if (verbose) {
    cat("=== Table Information: ", table_name, " ===\n")
  }

  # Get table count
  count_query <- paste0("SELECT COUNT(*) as count FROM ", table_name)
  record_count <- DBI::dbGetQuery(con, count_query)$count

  if (verbose) {
    cat("Total records:", record_count, "\n")
  }

  # Get schema information
  schema_info <- NULL
  if (include_schema) {
    tryCatch({
      schema_query <- paste0("PRAGMA table_info(", table_name, ")")
      schema_info <- DBI::dbGetQuery(con, schema_query)

      if (verbose && nrow(schema_info) > 0) {
        cat("\n=== Table Schema ===\n")
        cat(sprintf("%-20s %-15s %-10s %-10s\n", "Column", "Type", "Not Null", "Default"))
        cat(paste(rep("-", 60), collapse = ""), "\n")

        for (i in 1:nrow(schema_info)) {
          col <- schema_info[i, ]
          cat(sprintf("%-20s %-15s %-10s %-10s\n",
                      col$name,
                      col$type,
                      ifelse(col$notnull == 1, "YES", "NO"),
                      ifelse(is.na(col$dflt_value), "NULL", col$dflt_value)))
        }
      }
    }, error = function(e) {
      if (verbose) cat("Could not retrieve schema information\n")
    })
  }

  # Get sample data
  sample_data <- NULL
  if (record_count > 0 && sample_size > 0) {
    tryCatch({
      sample_query <- paste0("SELECT * FROM ", table_name, " LIMIT ", sample_size)
      sample_data <- DBI::dbGetQuery(con, sample_query)

      if (verbose) {
        cat("\n=== Sample Data ===\n")
        print(sample_data)
      }
    }, error = function(e) {
      if (verbose) cat("Could not retrieve sample data\n")
    })
  }

  # Get column statistics for numeric columns
  column_stats <- NULL
  if (record_count > 0 && !is.null(schema_info)) {
    numeric_columns <- schema_info$name[grepl("INTEGER|REAL|NUMERIC", schema_info$type, ignore.case = TRUE)]

    if (length(numeric_columns) > 0) {
      tryCatch({
        stats_queries <- lapply(numeric_columns, function(col) {
          paste0("SELECT '", col, "' as column_name, ",
                 "MIN(", col, ") as min_val, ",
                 "MAX(", col, ") as max_val, ",
                 "AVG(", col, ") as avg_val, ",
                 "COUNT(DISTINCT ", col, ") as unique_count ",
                 "FROM ", table_name, " WHERE ", col, " IS NOT NULL")
        })

        column_stats <- do.call(rbind, lapply(stats_queries, function(query) {
          DBI::dbGetQuery(con, query)
        }))

        if (verbose && nrow(column_stats) > 0) {
          cat("\n=== Numeric Column Statistics ===\n")
          cat(sprintf("%-15s %10s %10s %12s %10s\n", "Column", "Min", "Max", "Average", "Unique"))
          cat(paste(rep("-", 65), collapse = ""), "\n")

          for (i in 1:nrow(column_stats)) {
            stat <- column_stats[i, ]
            cat(sprintf("%-15s %10.2f %10.2f %12.2f %10d\n",
                        stat$column_name, stat$min_val, stat$max_val,
                        stat$avg_val, stat$unique_count))
          }
        }
      }, error = function(e) {
        if (verbose) cat("Could not calculate column statistics\n")
      })
    }
  }

  # Return comprehensive table info
  list(
    table_name = table_name,
    record_count = record_count,
    schema_info = schema_info,
    sample_data = sample_data,
    column_stats = column_stats,
    analysis_date = Sys.time()
  )
}

#' Quick table counts for all tables
#'
#' @param con A database connection object
#' @param verbose Logical. Print output. Default is TRUE
#'
#' @return Named vector of table counts
#' @export
get_table_counts <- function(con, verbose = TRUE) {

  if (!DBI::dbIsValid(con)) {
    stop("Invalid database connection")
  }

  all_tables <- DBI::dbListTables(con)
  counts <- setNames(numeric(length(all_tables)), all_tables)

  for (table in all_tables) {
    tryCatch({
      count_query <- paste0("SELECT COUNT(*) as count FROM ", table)
      counts[table] <- DBI::dbGetQuery(con, count_query)$count
    }, error = function(e) {
      counts[table] <- NA
    })
  }

  if (verbose) {
    cat("=== Table Record Counts ===\n")
    for (table in names(counts)) {
      cat(sprintf("%-25s: %s\n", table,
                  ifelse(is.na(counts[table]), "Error", format(counts[table], big.mark = ","))))
    }
  }

  return(counts)
}

#' New function to get detailed BLAST database information
#'
#' @param con A database connection object
#' @param db_name Optional. Filter by specific database name
#' @param verbose Logical. Print detailed output. Default is TRUE
#'
#' @return A list containing detailed database information
#' @export
get_blast_database_info <- function(con, db_name = NULL, verbose = TRUE) {
  if (!DBI::dbIsValid(con)) {
    stop("Invalid database connection")
  }

  # Check if metadata table exists
  tables <- DBI::dbListTables(con)
  if (!"blast_database_metadata" %in% tables) {
    if (verbose) cat("No BLAST database metadata available.\n")
    cat("Run perform_blast_db() with extract_db_metadata=TRUE to collect this information.\n")
    return(list())
  }

  # Build query
  base_query <- "
    SELECT bm.*, bp.project_id, p.project_name,
           COUNT(DISTINCT bp.blast_param_id) as usage_count,
           COUNT(DISTINCT br.blast_result_id) as total_results,
           COUNT(DISTINCT a.annotation_id) as total_annotations,
           MIN(bp.execution_date) as first_used,
           MAX(bp.execution_date) as last_used
    FROM blast_database_metadata bm
    JOIN blast_parameters bp ON bm.blast_param_id = bp.blast_param_id
    JOIN projects p ON bp.project_id = p.project_id
    LEFT JOIN blast_results br ON bp.blast_param_id = br.blast_param_id
    LEFT JOIN annotations a ON br.blast_result_id = a.blast_result_id
  "

  if (!is.null(db_name)) {
    query <- paste0(base_query, " WHERE bm.db_name = ? GROUP BY bm.metadata_id ORDER BY bm.extraction_date DESC")
    params <- list(db_name)
  } else {
    query <- paste0(base_query, " GROUP BY bm.metadata_id ORDER BY usage_count DESC, bm.extraction_date DESC")
    params <- list()
  }

  # Execute query
  db_info <- DBI::dbGetQuery(con, query, params = params)

  if (verbose && nrow(db_info) > 0) {
    cat("=== BLAST Database Information ===\n")

    for (i in 1:nrow(db_info)) {
      info <- db_info[i, ]
      cat(sprintf("\nDatabase: %s\n", info$db_name))
      cat(sprintf("Title: %s\n", info$db_title %||% "Unknown"))
      cat(sprintf("Path: %s\n", info$db_full_path))
      cat("Database Statistics:\n")
      cat(sprintf("  Sequences: %s\n", format(info$num_sequences %||% 0, big.mark = ",")))
      cat(sprintf("  Total length: %s bp\n", format(info$total_length %||% 0, big.mark = ",")))
      cat(sprintf("  Longest sequence: %s bp\n", format(info$longest_sequence %||% 0, big.mark = ",")))
      cat(sprintf("  Database date: %s\n", info$db_date %||% "Unknown"))
      if (!is.na(info$db_version) && !is.null(info$db_version)) {
        cat(sprintf("  Version: %s\n", info$db_version))
      }
      cat("Usage Statistics:\n")
      cat(sprintf("  Times used: %d\n", info$usage_count))
      cat(sprintf("  Total BLAST results: %s\n", format(info$total_results, big.mark = ",")))
      cat(sprintf("  Total annotations: %s\n", format(info$total_annotations, big.mark = ",")))
      cat(sprintf("  First used: %s\n", info$first_used))
      cat(sprintf("  Last used: %s\n", info$last_used))
      cat(sprintf("  Metadata extracted: %s\n", info$extraction_date))
      cat(paste(rep("-", 50), collapse = ""), "\n")
    }
  }

  return(db_info)
}
