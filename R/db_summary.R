#' Get comprehensive database summary
#'
#' This function provides a summary of the entire funseqR database contents
#'
#' @param con A database connection object
#' @param detail_level Character. Level of detail: "basic", "detailed", or "full". Default is "detailed"
#' @param include_samples Logical. Include sample data from tables. Default is FALSE
#' @param use_consolidated_chrom Logical. If TRUE, group minor scaffolds as "US" in chromosome distribution. Default is TRUE
#' @param verbose Logical. Print detailed output. Default is TRUE
#'
#' @return A list containing database summary information
#' @export
get_database_summary <- function(con, detail_level = "detailed", include_samples = FALSE, use_consolidated_chrom = TRUE, verbose = TRUE) {
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
  core_tables <- c("input_files", "vcf_data", "reference_genomes",
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

  # Get database-level summary (single project per database)
  database_overview <- list()
  if (verbose) cat("\n=== Database Overview ===\n")
  
  tryCatch({
    # Get file counts and analysis summary
    input_files_count <- table_counts[["input_files"]] %||% 0
    vcf_entries_count <- table_counts[["vcf_data"]] %||% 0
    blast_runs_count <- table_counts[["blast_parameters"]] %||% 0
    annotations_count <- table_counts[["annotations"]] %||% 0
    
    database_overview <- list(
      input_files = input_files_count,
      vcf_entries = vcf_entries_count,
      blast_runs = blast_runs_count,
      annotations = annotations_count
    )
    
    if (verbose) {
      cat(sprintf("Input files: %d\n", input_files_count))
      cat(sprintf("VCF entries: %d\n", vcf_entries_count))
      cat(sprintf("BLAST runs: %d\n", blast_runs_count))
      cat(sprintf("Annotations: %d\n", annotations_count))
    }
  }, error = function(e) {
    if (verbose) cat("Could not generate database overview\n")
  })

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

    # Chromosome distribution with optional consolidation
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
        
        # Apply consolidation if requested and main chromosomes are defined
        display_summary <- chrom_summary
        if (use_consolidated_chrom) {
          main_chroms <- get_main_chromosomes(con)
          if (!is.null(main_chroms) && length(main_chroms) > 0) {
            # Create consolidated data for display
            vcf_data_temp <- data.frame(chromosome = rep(chrom_summary$chromosome, chrom_summary$variant_count))
            consolidated_temp <- consolidate_scaffolds(vcf_data_temp, main_chroms, verbose = FALSE)
            display_summary <- aggregate(variant_count ~ chromosome, 
                                       data = data.frame(chromosome = consolidated_temp$chromosome, 
                                                        variant_count = 1), 
                                       FUN = sum)
            display_summary <- display_summary[order(display_summary$variant_count, decreasing = TRUE), ]
          }
        }
        
        detailed_info$chromosome_distribution <- chrom_summary  # Store original
        detailed_info$chromosome_distribution_display <- display_summary  # Store display version

        if (verbose && nrow(display_summary) > 0) {
          cat("\n=== Variants by Chromosome ===\n")
          if (use_consolidated_chrom) {
            main_chroms <- get_main_chromosomes(con)
            if (!is.null(main_chroms) && length(main_chroms) > 0) {
              cat("(Using consolidated view - minor scaffolds grouped as 'US')\n")
            }
          }
          
          total_variants <- sum(display_summary$variant_count)
          display_limit <- min(nrow(display_summary), 15)  # Show more for consolidated view
          
          for (i in 1:display_limit) {
            cs <- display_summary[i, ]
            pct <- round(100 * cs$variant_count / total_variants, 1)
            cat(sprintf("%-10s: %s variants (%s%%)\n", cs$chromosome, 
                       format(cs$variant_count, big.mark = ","), pct))
          }
          if (nrow(display_summary) > display_limit) {
            cat(sprintf("... and %d more chromosomes\n", nrow(display_summary) - display_limit))
          }
          
          if (use_consolidated_chrom && !is.null(get_main_chromosomes(con))) {
            cat(sprintf("\nOriginal chromosome count: %d, Consolidated view: %d\n", 
                       nrow(chrom_summary), nrow(display_summary)))
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
    sample_tables <- c("input_files", "vcf_data", "annotations", "go_terms")

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
    database_overview = database_overview,
    detailed_info = detailed_info,
    sample_data = sample_data,
    summary_date = Sys.time()
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
                 "COUNT(DISTINCT ", col, ") as unique_count ",
                 "FROM ", table_name, " WHERE ", col, " IS NOT NULL")
        })

        column_stats <- do.call(rbind, lapply(stats_queries, function(query) {
          DBI::dbGetQuery(con, query)
        }))

        if (verbose && nrow(column_stats) > 0) {
          cat("\n=== Numeric Column Statistics ===\n")
          cat(sprintf("%-15s %10s %10s %10s\n", "Column", "Min", "Max", "Unique"))
          cat(paste(rep("-", 50), collapse = ""), "\n")

          for (i in 1:nrow(column_stats)) {
            stat <- column_stats[i, ]
            cat(sprintf("%-15s %10.2f %10.2f %10d\n",
                        stat$column_name, stat$min_val, stat$max_val,
                        stat$unique_count))
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

  # Build query (simplified for single-project design)
  base_query <- "
    SELECT bm.*,
           COUNT(DISTINCT bp.blast_param_id) as usage_count,
           COUNT(DISTINCT br.blast_result_id) as total_results,
           COUNT(DISTINCT a.annotation_id) as total_annotations,
           MIN(bp.execution_date) as first_used,
           MAX(bp.execution_date) as last_used
    FROM blast_database_metadata bm
    JOIN blast_parameters bp ON bm.blast_param_id = bp.blast_param_id
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
