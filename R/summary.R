# EXPORTED

#' Database summary functions for funseqR
#'
#' These functions provide summary information about the contents of a funseqR database.
#'

#' Summarize funseqR database contents
#'
#' This function provides comprehensive summaries of funseqR database contents.
#' Supports multiple summary types organized by workflow stages.
#'
#' @param con A database connection object created by \code{\link{connect_funseq_db}} or \code{\link{create_funseq_db}}.
#' @param type Character string specifying the type of summary to generate:
#'   \describe{
#'     \item{database}{Overall database summary with record counts for all tables}
#'     \item{input_data}{Detailed summaries for core data tables (input_files, vcf_data, reference_genomes, reference_sequences, flanking_sequences)}
#'     \item{blast}{BLAST-related summaries (blast_parameters, blast_results, blast_database_metadata)}
#'     \item{annotation}{Annotation summaries (annotations, go_terms, kegg_references, pfam_domains, interpro_families, eggnog_categories)}
#'     \item{analyses}{Analysis summaries (ora_analyses, ora_results, candidate_loci, locus_statistics)}
#'     \item{methods}{Methods documentation for reproducibility (logged commands, parameters, software versions)}
#'   }
#'
#' @return For type "database": data.frame with table_name and record_count columns.
#'   For other types: named list with detailed data.frame summaries for each table in that group.
#'
#' @details
#' The function provides summaries organized by workflow stages:
#' \itemize{
#'   \item \strong{database}: Template-based overview of all tables with record counts
#'   \item \strong{input_data}: File information, VCF statistics, genome details, sequence info
#'   \item \strong{blast}: BLAST run parameters, hit statistics, database metadata
#'   \item \strong{annotation}: Functional annotation coverage and distributions
#'   \item \strong{analyses}: Enrichment results and candidate loci information
#'   \item \strong{methods}: Complete methods documentation for reproducibility including commands, parameters, and software versions
#' }
#'
#' @importFrom DBI dbListTables dbGetQuery
#'
#' @examples
#' \dontrun{
#' # Connect to a funseqR database
#' con <- connect_funseq_db("my_project.db")
#' 
#' # Get overall database summary
#' db_summary <- funseqR_summary(con, type = "database")
#' print(db_summary)
#' 
#' # Get detailed input data summaries
#' input_summary <- funseqR_summary(con, type = "input_data")
#' names(input_summary)  # "input_files", "vcf_data", etc.
#' 
#' # Get BLAST-related summaries
#' blast_summary <- funseqR_summary(con, type = "blast")
#' 
#' # Get annotation summaries
#' annotation_summary <- funseqR_summary(con, type = "annotation")
#' 
#' # Get analysis summaries
#' analysis_summary <- funseqR_summary(con, type = "analyses")
#' 
#' # Get methods documentation for reproducibility
#' methods_summary <- funseqR_summary(con, type = "methods")
#' 
#' # Close connection
#' close_funseq_db(con)
#' }
#'
#' @export
funseqR_summary <- function(con, type = "database") {
  # Validate connection
  if (!DBI::dbIsValid(con)) {
    stop("Invalid database connection.")
  }
  
  # Validate type parameter
  valid_types <- c("database", "input_data", "blast", "annotation", "analyses", "methods")
  if (!type %in% valid_types) {
    stop("Invalid type. Supported types: ", paste(valid_types, collapse = ", "))
  }
  
  if (type == "database") {
    # Pre-defined template with all possible funseqR tables
    # Listed in logical order: input -> processing -> analysis -> output
    template <- data.frame(
      table_name = c(
        # Core data tables
        "input_files",
        "vcf_data", 
        "reference_genomes",
        "reference_sequences",
        "flanking_sequences",
        
        # BLAST-related tables
        "blast_parameters",
        "blast_results",
        "blast_database_metadata",
        
        # Annotation tables
        "annotations",
        "go_terms",
        "kegg_references", 
        "pfam_domains",
        "interpro_families",
        "eggnog_categories",
        
        # Analysis tables
        "ora_analyses",
        "ora_results",
        "analysis_reports",
        
        # Utility tables
        "uniprot_cache",
        "locus_statistics",
        "candidate_loci",
        "method_log"
      ),
      record_count = 0L,
      stringsAsFactors = FALSE
    )
    
    # Get list of existing tables
    existing_tables <- DBI::dbListTables(con)
    
    # Count records for each existing table
    for (i in seq_len(nrow(template))) {
      table_name <- template$table_name[i]
      if (table_name %in% existing_tables) {
        tryCatch({
          count_result <- DBI::dbGetQuery(con, 
            paste("SELECT COUNT(*) as count FROM", table_name))
          template$record_count[i] <- count_result$count
        }, error = function(e) {
          # If there's an error querying the table, leave count as 0
          warning("Could not query table '", table_name, "': ", e$message)
        })
      }
    }
    
    return(template)
    
  } else if (type == "input_data") {
    return(.get_input_data_summary(con))
    
  } else if (type == "blast") {
    return(.get_blast_summary(con))
    
  } else if (type == "annotation") {
    return(.get_annotation_summary(con))
    
  } else if (type == "analyses") {
    return(.get_analyses_summary(con))
    
  } else if (type == "methods") {
    return(.get_methods_summary(con))
  }
}

# ==========================
# HELPER FUNCTIONS FOR GROUP SUMMARIES
# ==========================

#' Get input data summaries
#' @keywords internal
.get_input_data_summary <- function(con) {
  existing_tables <- DBI::dbListTables(con)
  result <- list()
  
  # input_files summary
  if ("input_files" %in% existing_tables) {
    result$input_files <- tryCatch({
      DBI::dbGetQuery(con, "
        SELECT 
          file_id,
          file_type,
          file_name,
          ROUND(LENGTH(file_hash) / 2.0 / 1024 / 1024, 2) as estimated_size_mb,
          import_date
        FROM input_files 
        ORDER BY file_id
      ")
    }, error = function(e) {
      data.frame(message = paste("Error querying input_files:", e$message))
    })
  } else {
    result$input_files <- data.frame(message = "Table not found")
  }
  
  # vcf_data summary
  if ("vcf_data" %in% existing_tables) {
    result$vcf_data <- tryCatch({
      DBI::dbGetQuery(con, "
        SELECT 
          chromosome,
          COUNT(*) as variant_count,
          MIN(position) as min_position,
          MAX(position) as max_position,
          ROUND(AVG(qual), 2) as avg_quality
        FROM vcf_data 
        GROUP BY chromosome 
        ORDER BY chromosome
      ")
    }, error = function(e) {
      data.frame(message = paste("Error querying vcf_data:", e$message))
    })
  } else {
    result$vcf_data <- data.frame(message = "Table not found")
  }
  
  # reference_genomes summary
  if ("reference_genomes" %in% existing_tables) {
    result$reference_genomes <- tryCatch({
      DBI::dbGetQuery(con, "
        SELECT 
          rg.genome_id,
          rg.genome_name,
          rg.genome_build,
          COUNT(rs.sequence_id) as sequence_count,
          ROUND(SUM(rs.sequence_length) / 1000000.0, 2) as total_length_mb
        FROM reference_genomes rg
        LEFT JOIN reference_sequences rs ON rg.genome_id = rs.genome_id
        GROUP BY rg.genome_id, rg.genome_name, rg.genome_build
        ORDER BY rg.genome_id
      ")
    }, error = function(e) {
      data.frame(message = paste("Error querying reference_genomes:", e$message))
    })
  } else {
    result$reference_genomes <- data.frame(message = "Table not found")
  }
  
  # reference_sequences summary
  if ("reference_sequences" %in% existing_tables) {
    result$reference_sequences <- tryCatch({
      DBI::dbGetQuery(con, "
        SELECT 
          'All sequences' as sequence_type,
          COUNT(*) as count,
          ROUND(AVG(sequence_length), 0) as avg_length,
          MAX(sequence_length) as max_length,
          MIN(sequence_length) as min_length
        FROM reference_sequences
        UNION ALL
        SELECT 
          CASE 
            WHEN sequence_length > 1000000 THEN 'Large (>1Mb)'
            WHEN sequence_length > 100000 THEN 'Medium (100Kb-1Mb)'
            ELSE 'Small (<100Kb)'
          END as sequence_type,
          COUNT(*) as count,
          ROUND(AVG(sequence_length), 0) as avg_length,
          MAX(sequence_length) as max_length,
          MIN(sequence_length) as min_length
        FROM reference_sequences
        GROUP BY 
          CASE 
            WHEN sequence_length > 1000000 THEN 'Large (>1Mb)'
            WHEN sequence_length > 100000 THEN 'Medium (100Kb-1Mb)'
            ELSE 'Small (<100Kb)'
          END
      ")
    }, error = function(e) {
      data.frame(message = paste("Error querying reference_sequences:", e$message))
    })
  } else {
    result$reference_sequences <- data.frame(message = "Table not found")
  }
  
  # flanking_sequences summary  
  if ("flanking_sequences" %in% existing_tables) {
    result$flanking_sequences <- tryCatch({
      DBI::dbGetQuery(con, "
        SELECT 
          seq_type,
          COUNT(*) as count,
          flank_size,
          ROUND(AVG(seq_length), 0) as avg_seq_length,
          ROUND(AVG(end_position - start_position), 0) as avg_extracted_length
        FROM flanking_sequences
        GROUP BY seq_type, flank_size
        ORDER BY seq_type, flank_size
      ")
    }, error = function(e) {
      data.frame(message = paste("Error querying flanking_sequences:", e$message))
    })
  } else {
    result$flanking_sequences <- data.frame(message = "Table not found")
  }
  
  return(result)
}

#' Get BLAST summaries
#' @keywords internal
.get_blast_summary <- function(con) {
  existing_tables <- DBI::dbListTables(con)
  result <- list()
  
  # blast_parameters summary
  if ("blast_parameters" %in% existing_tables) {
    result$blast_parameters <- tryCatch({
      DBI::dbGetQuery(con, "
        SELECT 
          blast_param_id,
          blast_type,
          db_name,
          e_value,
          max_hits,
          execution_date
        FROM blast_parameters 
        ORDER BY blast_param_id
      ")
    }, error = function(e) {
      data.frame(message = paste("Error querying blast_parameters:", e$message))
    })
  } else {
    result$blast_parameters <- data.frame(message = "Table not found")
  }
  
  # blast_results summary
  if ("blast_results" %in% existing_tables) {
    result$blast_results <- tryCatch({
      DBI::dbGetQuery(con, "
        SELECT 
          bp.blast_type,
          COUNT(*) as total_hits,
          ROUND(AVG(br.percent_identity), 2) as avg_identity,
          ROUND(MIN(br.percent_identity), 2) as min_identity,
          ROUND(MAX(br.percent_identity), 2) as max_identity,
          ROUND(AVG(br.e_value), 6) as avg_e_value,
          COUNT(DISTINCT br.hit_accession) as unique_hits
        FROM blast_results br
        JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
        GROUP BY bp.blast_type
        ORDER BY bp.blast_type
      ")
    }, error = function(e) {
      data.frame(message = paste("Error querying blast_results:", e$message))
    })
  } else {
    result$blast_results <- data.frame(message = "Table not found")
  }
  
  # blast_database_metadata summary
  if ("blast_database_metadata" %in% existing_tables) {
    result$blast_database_metadata <- tryCatch({
      DBI::dbGetQuery(con, "
        SELECT 
          db_name,
          db_title,
          num_sequences,
          ROUND(total_length / 1000000.0, 2) as total_length_mb,
          db_date,
          db_version,
          extraction_date
        FROM blast_database_metadata 
        ORDER BY db_name
      ")
    }, error = function(e) {
      data.frame(message = paste("Error querying blast_database_metadata:", e$message))
    })
  } else {
    result$blast_database_metadata <- data.frame(message = "Table not found")
  }
  
  return(result)
}

#' Get annotation summaries
#' @keywords internal
.get_annotation_summary <- function(con) {
  existing_tables <- DBI::dbListTables(con)
  result <- list()
  
  # annotations summary
  if ("annotations" %in% existing_tables) {
    result$annotations <- tryCatch({
      DBI::dbGetQuery(con, "
        SELECT 
          COUNT(*) as total_annotations,
          COUNT(DISTINCT uniprot_accession) as unique_proteins,
          COUNT(DISTINCT br.flanking_id) as annotated_loci,
          MIN(a.retrieval_date) as first_retrieval,
          MAX(a.retrieval_date) as last_retrieval
        FROM annotations a
        JOIN blast_results br ON a.blast_result_id = br.blast_result_id
      ")
    }, error = function(e) {
      data.frame(message = paste("Error querying annotations:", e$message))
    })
  } else {
    result$annotations <- data.frame(message = "Table not found")
  }
  
  # go_terms summary
  if ("go_terms" %in% existing_tables) {
    result$go_terms <- tryCatch({
      DBI::dbGetQuery(con, "
        SELECT 
          go_category,
          COUNT(*) as term_count,
          COUNT(DISTINCT go_id) as unique_terms,
          COUNT(DISTINCT annotation_id) as annotated_proteins
        FROM go_terms
        GROUP BY go_category
        ORDER BY go_category
      ")
    }, error = function(e) {
      data.frame(message = paste("Error querying go_terms:", e$message))
    })
  } else {
    result$go_terms <- data.frame(message = "Table not found")
  }
  
  # kegg_references summary
  if ("kegg_references" %in% existing_tables) {
    result$kegg_references <- tryCatch({
      DBI::dbGetQuery(con, "
        SELECT 
          'KEGG Pathways' as category,
          COUNT(*) as total_associations,
          COUNT(DISTINCT kegg_id) as unique_pathways,
          COUNT(DISTINCT annotation_id) as annotated_proteins,
          COUNT(DISTINCT CASE WHEN pathway_name IS NOT NULL AND pathway_name != '' THEN kegg_id END) as pathways_with_names
        FROM kegg_references
      ")
    }, error = function(e) {
      data.frame(message = paste("Error querying kegg_references:", e$message))
    })
  } else {
    result$kegg_references <- data.frame(message = "Table not found")
  }
  
  # pfam_domains summary
  if ("pfam_domains" %in% existing_tables) {
    result$pfam_domains <- tryCatch({
      DBI::dbGetQuery(con, "
        SELECT 
          match_status,
          COUNT(*) as association_count,
          COUNT(DISTINCT pfam_id) as unique_domains,
          COUNT(DISTINCT annotation_id) as annotated_proteins
        FROM pfam_domains
        GROUP BY match_status
        ORDER BY match_status
      ")
    }, error = function(e) {
      data.frame(message = paste("Error querying pfam_domains:", e$message))
    })
  } else {
    result$pfam_domains <- data.frame(message = "Table not found")
  }
  
  # interpro_families summary
  if ("interpro_families" %in% existing_tables) {
    result$interpro_families <- tryCatch({
      DBI::dbGetQuery(con, "
        SELECT 
          'InterPro Families' as category,
          COUNT(*) as total_associations,
          COUNT(DISTINCT interpro_id) as unique_families,
          COUNT(DISTINCT annotation_id) as annotated_proteins,
          COUNT(DISTINCT CASE WHEN family_name IS NOT NULL AND family_name != '' THEN interpro_id END) as families_with_names
        FROM interpro_families
      ")
    }, error = function(e) {
      data.frame(message = paste("Error querying interpro_families:", e$message))
    })
  } else {
    result$interpro_families <- data.frame(message = "Table not found")
  }
  
  # eggnog_categories summary
  if ("eggnog_categories" %in% existing_tables) {
    result$eggnog_categories <- tryCatch({
      DBI::dbGetQuery(con, "
        SELECT 
          taxonomic_scope,
          COUNT(*) as association_count,
          COUNT(DISTINCT eggnog_id) as unique_categories,
          COUNT(DISTINCT annotation_id) as annotated_proteins
        FROM eggnog_categories
        GROUP BY taxonomic_scope
        ORDER BY taxonomic_scope
      ")
    }, error = function(e) {
      data.frame(message = paste("Error querying eggnog_categories:", e$message))
    })
  } else {
    result$eggnog_categories <- data.frame(message = "Table not found")
  }
  
  return(result)
}

#' Get analyses summaries
#' @keywords internal
.get_analyses_summary <- function(con) {
  existing_tables <- DBI::dbListTables(con)
  result <- list()
  
  # ora_analyses summary
  if ("ora_analyses" %in% existing_tables) {
    result$ora_analyses <- tryCatch({
      DBI::dbGetQuery(con, "
        SELECT 
          analysis_id,
          annotation_type,
          term_type,
          analysis_date,
          total_foreground_genes,
          total_background_genes,
          enrichment_method
        FROM ora_analyses 
        ORDER BY analysis_id
      ")
    }, error = function(e) {
      data.frame(message = paste("Error querying ora_analyses:", e$message))
    })
  } else {
    result$ora_analyses <- data.frame(message = "Table not found")
  }
  
  # ora_results summary
  if ("ora_results" %in% existing_tables) {
    result$ora_results <- tryCatch({
      DBI::dbGetQuery(con, "
        SELECT 
          oa.annotation_type,
          oa.term_type,
          COUNT(*) as total_terms_tested,
          COUNT(CASE WHEN ore.p_adjusted <= 0.05 THEN 1 END) as significant_05,
          COUNT(CASE WHEN ore.p_adjusted <= 0.1 THEN 1 END) as significant_10,
          ROUND(MIN(ore.p_adjusted), 6) as min_fdr,
          ROUND(MAX(ore.fold_enrichment), 2) as max_enrichment,
          ROUND(AVG(ore.fold_enrichment), 2) as avg_enrichment
        FROM ora_results ore
        JOIN ora_analyses oa ON ore.analysis_id = oa.analysis_id
        GROUP BY oa.annotation_type, oa.term_type
        ORDER BY oa.annotation_type, oa.term_type
      ")
    }, error = function(e) {
      data.frame(message = paste("Error querying ora_results:", e$message))
    })
  } else {
    result$ora_results <- data.frame(message = "Table not found")
  }
  
  # candidate_loci summary
  if ("candidate_loci" %in% existing_tables) {
    result$candidate_loci <- tryCatch({
      DBI::dbGetQuery(con, "
        SELECT 
          method,
          COUNT(*) as loci_count,
          ROUND(AVG(threshold), 4) as avg_threshold,
          MIN(threshold) as min_threshold,
          MAX(threshold) as max_threshold,
          COUNT(DISTINCT chromosome) as chromosomes_with_candidates
        FROM candidate_loci
        GROUP BY method
        ORDER BY method
      ")
    }, error = function(e) {
      data.frame(message = paste("Error querying candidate_loci:", e$message))
    })
  } else {
    result$candidate_loci <- data.frame(message = "Table not found")
  }
  
  # locus_statistics summary
  if ("locus_statistics" %in% existing_tables) {
    result$locus_statistics <- tryCatch({
      DBI::dbGetQuery(con, "
        SELECT 
          COUNT(*) as total_loci,
          COUNT(DISTINCT chromosome) as chromosomes,
          ROUND(MIN(statistic), 4) as min_statistic,
          ROUND(MAX(statistic), 4) as max_statistic,
          ROUND(AVG(statistic), 4) as mean_statistic,
          ROUND(
            (SELECT statistic FROM locus_statistics ORDER BY statistic LIMIT 1 OFFSET (
              (SELECT COUNT(*) FROM locus_statistics)/2
            ))
          , 4) as median_statistic
        FROM locus_statistics
      ")
    }, error = function(e) {
      data.frame(message = paste("Error querying locus_statistics:", e$message))
    })
  } else {
    result$locus_statistics <- data.frame(message = "Table not found")
  }
  
  return(result)
}

#' Get methods summaries for reproducibility
#' @keywords internal
.get_methods_summary <- function(con) {
  existing_tables <- DBI::dbListTables(con)
  result <- list()
  
  # Check if method_log table exists
  if (!"method_log" %in% existing_tables) {
    result$status <- data.frame(
      message = "Method logging not available - method_log table not found",
      note = "This database was created before method logging was implemented"
    )
    return(result)
  }
  
  # Get overall method log summary
  result$overview <- tryCatch({
    DBI::dbGetQuery(con, "
      SELECT 
        method_type,
        COUNT(*) as total_executions,
        COUNT(DISTINCT function_name) as unique_functions,
        COUNT(CASE WHEN success = 1 THEN 1 END) as successful_executions,
        COUNT(CASE WHEN success = 0 THEN 1 END) as failed_executions,
        MIN(execution_date) as first_execution,
        MAX(execution_date) as last_execution,
        ROUND(AVG(execution_time_seconds), 2) as avg_execution_time_seconds
      FROM method_log
      GROUP BY method_type
      ORDER BY method_type
    ")
  }, error = function(e) {
    data.frame(message = paste("Error querying method log overview:", e$message))
  })
  
  # Get BLAST methods
  result$blast_methods <- tryCatch({
    DBI::dbGetQuery(con, "
      SELECT 
        function_name,
        command_text,
        parameters_json,
        execution_date,
        execution_time_seconds,
        success
      FROM method_log
      WHERE method_type = 'blast'
      ORDER BY execution_date DESC
    ")
  }, error = function(e) {
    data.frame(message = paste("Error querying BLAST methods:", e$message))
  })
  
  # Get annotation methods
  result$annotation_methods <- tryCatch({
    DBI::dbGetQuery(con, "
      SELECT 
        function_name,
        command_text,
        parameters_json,
        execution_date,
        execution_time_seconds,
        success
      FROM method_log
      WHERE method_type = 'annotation'
      ORDER BY execution_date DESC
    ")
  }, error = function(e) {
    data.frame(message = paste("Error querying annotation methods:", e$message))
  })
  
  # Get enrichment methods
  result$enrichment_methods <- tryCatch({
    DBI::dbGetQuery(con, "
      SELECT 
        function_name,
        command_text,
        parameters_json,
        execution_date,
        execution_time_seconds,
        success
      FROM method_log
      WHERE method_type = 'enrichment'
      ORDER BY execution_date DESC
    ")
  }, error = function(e) {
    data.frame(message = paste("Error querying enrichment methods:", e$message))
  })
  
  # Get sequence processing methods
  result$sequence_methods <- tryCatch({
    DBI::dbGetQuery(con, "
      SELECT 
        function_name,
        command_text,
        parameters_json,
        execution_date,
        execution_time_seconds,
        success
      FROM method_log
      WHERE method_type = 'sequence'
      ORDER BY execution_date DESC
    ")
  }, error = function(e) {
    data.frame(message = paste("Error querying sequence methods:", e$message))
  })
  
  # Get statistical methods
  result$statistical_methods <- tryCatch({
    DBI::dbGetQuery(con, "
      SELECT 
        function_name,
        command_text,
        parameters_json,
        execution_date,
        execution_time_seconds,
        success
      FROM method_log
      WHERE method_type = 'statistical'
      ORDER BY execution_date DESC
    ")
  }, error = function(e) {
    data.frame(message = paste("Error querying statistical methods:", e$message))
  })
  
  # Get software environment information
  result$software_environment <- tryCatch({
    DBI::dbGetQuery(con, "
      SELECT DISTINCT
        r_version,
        package_versions,
        MIN(execution_date) as first_used,
        MAX(execution_date) as last_used
      FROM method_log
      WHERE r_version IS NOT NULL
      GROUP BY r_version, package_versions
      ORDER BY first_used DESC
    ")
  }, error = function(e) {
    data.frame(message = paste("Error querying software environment:", e$message))
  })
  
  # Generate formatted methods text for scientific papers
  result$methods_text <- tryCatch({
    .format_methods_text(con)
  }, error = function(e) {
    data.frame(message = paste("Error generating methods text:", e$message))
  })
  
  return(result)
}

#' Format methods information for scientific papers
#' @keywords internal
.format_methods_text <- function(con) {
  # Get summary statistics for formatting
  overview <- DBI::dbGetQuery(con, "
    SELECT 
      method_type,
      COUNT(*) as executions,
      COUNT(DISTINCT function_name) as functions
    FROM method_log
    GROUP BY method_type
  ")
  
  # Get software versions
  software <- DBI::dbGetQuery(con, "
    SELECT DISTINCT r_version, package_versions
    FROM method_log
    WHERE r_version IS NOT NULL
    LIMIT 1
  ")
  
  # Get BLAST parameters
  blast_params <- DBI::dbGetQuery(con, "
    SELECT DISTINCT parameters_json
    FROM method_log
    WHERE method_type = 'blast'
    LIMIT 1
  ")
  
  # Get annotation parameters
  annotation_params <- DBI::dbGetQuery(con, "
    SELECT DISTINCT parameters_json
    FROM method_log
    WHERE method_type = 'annotation'
    LIMIT 1
  ")
  
  # Get enrichment parameters
  enrichment_params <- DBI::dbGetQuery(con, "
    SELECT DISTINCT parameters_json
    FROM method_log
    WHERE method_type = 'enrichment'
    LIMIT 1
  ")
  
  # Format the text
  methods_text <- ""
  
  # Software section
  if (nrow(software) > 0 && !is.na(software$r_version[1])) {
    methods_text <- paste0(methods_text, 
      "All analyses were performed using the funseqR package in R version ", 
      software$r_version[1], ". ")
  }
  
  # BLAST section
  if ("blast" %in% overview$method_type) {
    methods_text <- paste0(methods_text, 
      "Sequence similarity searches were performed using BLAST. ")
    if (nrow(blast_params) > 0 && !is.na(blast_params$parameters_json[1])) {
      methods_text <- paste0(methods_text, "BLAST parameters: ", 
                            blast_params$parameters_json[1], ". ")
    }
  }
  
  # Annotation section
  if ("annotation" %in% overview$method_type) {
    methods_text <- paste0(methods_text, 
      "Functional annotations were retrieved from the UniProt database. ")
    if (nrow(annotation_params) > 0 && !is.na(annotation_params$parameters_json[1])) {
      methods_text <- paste0(methods_text, "Annotation parameters: ", 
                            annotation_params$parameters_json[1], ". ")
    }
  }
  
  # Enrichment section
  if ("enrichment" %in% overview$method_type) {
    methods_text <- paste0(methods_text, 
      "Over-representation analysis was performed using the clusterProfiler package. ")
    if (nrow(enrichment_params) > 0 && !is.na(enrichment_params$parameters_json[1])) {
      methods_text <- paste0(methods_text, "Enrichment parameters: ", 
                            enrichment_params$parameters_json[1], ". ")
    }
  }
  
  return(data.frame(
    section = "Methods",
    text = methods_text,
    note = "This text was automatically generated from logged method parameters"
  ))
}