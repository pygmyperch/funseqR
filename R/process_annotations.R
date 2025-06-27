#' Process Functional Annotations into Standardized Summary Table
#'
#' Processes functional annotation results from annotate_blast_results() into a
#' comprehensive, standardized data frame linking each locus to functional annotations.
#' Automatically applies name standardization using KEGGREST and includes KEGG BRITE
#' hierarchy and module classifications.
#'
#' @param con Database connection object from annotate_blast_results()
#' @param include Character vector. Types of annotations to include: "GO", "KEGG", or both. Default is c("GO", "KEGG")
#' @param blast_param_id Integer. Optional. Specific BLAST parameter set to use. If NULL, uses all available annotations
#' @param export_csv Character. Optional. File path to export results as CSV
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return Data frame with standardized functional annotations:
#' \itemize{
#'   \item locus_id: Unique identifier for each genomic locus
#'   \item chromosome: Chromosome name
#'   \item position: Genomic position
#'   \item gene_name: Primary gene name (from UniProt)
#'   \item protein_name: Protein description
#'   \item uniprot_accession: UniProt accession number
#'   \item go_terms: Semi-colon separated GO term IDs (if "GO" included)
#'   \item go_names: Semi-colon separated GO term names (if "GO" included)
#'   \item go_categories: Semi-colon separated GO categories (BP/MF/CC) (if "GO" included)
#'   \item kegg_pathways: Semi-colon separated KEGG pathway IDs (if "KEGG" included)
#'   \item kegg_pathway_names: Semi-colon separated KEGG pathway names (if "KEGG" included)
#'   \item kegg_brite_categories: Semi-colon separated KEGG BRITE classifications (if "KEGG" included)
#'   \item kegg_modules: Semi-colon separated KEGG module assignments (if "KEGG" included)
#' }
#'
#' @details
#' This function streamlines the functional annotation workflow by:
#' 
#' \enumerate{
#'   \item \strong{Data Integration}: Combines genomic positions with functional annotations
#'   \item \strong{GO Processing}: Extracts GO terms, names, and ontology categories
#'   \item \strong{KEGG Enhancement}: Uses KEGGREST API to improve pathway names and classifications
#'   \item \strong{KEGG BRITE Classification}: Maps pathways to official KEGG BRITE hierarchy
#'   \item \strong{Module Assignment}: Assigns pathways to functional modules
#'   \item \strong{Standardized Output}: Creates one comprehensive table for all downstream analyses
#' }
#' 
#' The resulting data frame can be used directly with enrichment analysis functions
#' or exported for external analysis tools.
#'
#' @examples
#' \dontrun{
#' # Complete workflow
#' blast_results <- run_blast(sequences, database)
#' con <- annotate_blast_results(blast_results)
#' 
#' # Process all annotations
#' annotations <- process_annotations(con, include = c("GO", "KEGG"))
#' 
#' # Process only KEGG annotations
#' kegg_only <- process_annotations(con, include = "KEGG")
#' 
#' # Export to CSV for external analysis
#' annotations <- process_annotations(
#'   con, 
#'   include = c("GO", "KEGG"),
#'   export_csv = "functional_annotations.csv"
#' )
#' 
#' # Use with enrichment analysis
#' go_enrichment <- run_go_enrichment_analysis(annotations, candidate_loci)
#' kegg_enrichment <- run_kegg_enrichment_analysis(annotations, candidate_loci)
#' }
#'
#' @importFrom dplyr group_by summarise first rowwise ungroup
#' @importFrom magrittr %>%
#' @export
process_annotations <- function(con, include = c("GO", "KEGG", "Pfam", "InterPro", "eggNOG"), blast_param_id = NULL, 
                               export_csv = NULL, verbose = TRUE) {
  
  if (verbose) message("Processing functional annotations...")
  
  # Validate inputs
  valid_types <- c("GO", "KEGG", "Pfam", "InterPro", "eggNOG")
  include <- match.arg(include, valid_types, several.ok = TRUE)
  
  if (verbose) {
    message("  - Including annotations: ", paste(include, collapse = ", "))
    if (!is.null(blast_param_id)) {
      message("  - Using BLAST parameter set: ", blast_param_id)
    }
  }
  
  # Check database tables
  tables <- DBI::dbListTables(con)
  required_base_tables <- c("vcf_data", "flanking_sequences", "blast_results", "annotations")
  missing_tables <- required_base_tables[!required_base_tables %in% tables]
  
  if (length(missing_tables) > 0) {
    stop("Missing required tables: ", paste(missing_tables, collapse = ", "))
  }
  
  # Build base query for loci with annotations
  base_conditions <- c()
  params <- list()
  
  if (!is.null(blast_param_id)) {
    base_conditions <- c(base_conditions, "bp.blast_param_id = ?")
    params <- c(params, list(blast_param_id))
  }
  
  base_where <- if (length(base_conditions) > 0) {
    paste("WHERE", paste(base_conditions, collapse = " AND "))
  } else {
    ""
  }
  
  # Start with base loci and annotation information
  if (verbose) message("  - Extracting basic annotation data...")
  
  base_query <- paste0("
    SELECT DISTINCT
      vd.vcf_id || '_' || vd.chromosome || '_' || vd.position as locus_id,
      vd.chromosome,
      vd.position,
      a.annotation_id,
      a.uniprot_accession,
      a.entry_name,
      a.gene_names,
      br.blast_result_id
    FROM vcf_data vd
    JOIN flanking_sequences fs ON vd.vcf_id = fs.vcf_id
    JOIN blast_results br ON fs.flanking_id = br.flanking_id
    JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
    JOIN annotations a ON br.blast_result_id = a.blast_result_id
    ", base_where, "
    ORDER BY vd.chromosome, vd.position
  ")
  
  base_data <- if (length(params) > 0) {
    DBI::dbGetQuery(con, base_query, params)
  } else {
    DBI::dbGetQuery(con, base_query)
  }
  
  if (nrow(base_data) == 0) {
    warning("No annotation data found with current filters")
    return(data.frame())
  }
  
  if (verbose) message("    - Found ", nrow(base_data), " annotated loci")
  
  # Process basic information
  result_data <- base_data %>%
    group_by(locus_id, chromosome, position) %>%
    summarise(
      gene_name = first(na.omit(gene_names))[1],
      protein_name = first(na.omit(entry_name))[1],  # Use entry_name as protein description
      uniprot_accession = paste(unique(uniprot_accession), collapse = ";"),
      annotation_ids = list(unique(annotation_id)),
      .groups = "drop"
    )
  
  # Process GO annotations if requested
  if ("GO" %in% include && "go_terms" %in% tables) {
    if (verbose) message("  - Processing GO annotations...")
    result_data <- .process_go_annotations(con, result_data, base_where, params, verbose)
  }
  
  # Process KEGG annotations if requested  
  if ("KEGG" %in% include && "kegg_references" %in% tables) {
    if (verbose) message("  - Processing KEGG annotations with enhancements...")
    result_data <- .process_kegg_annotations(con, result_data, base_where, params, verbose)
  }
  
  # Process Pfam annotations if requested
  if ("Pfam" %in% include && "pfam_domains" %in% tables) {
    if (verbose) message("  - Processing Pfam domain annotations...")
    result_data <- .process_pfam_annotations(con, result_data, base_where, params, verbose)
  }
  
  # Process InterPro annotations if requested
  if ("InterPro" %in% include && "interpro_families" %in% tables) {
    if (verbose) message("  - Processing InterPro family annotations...")
    result_data <- .process_interpro_annotations(con, result_data, base_where, params, verbose)
  }
  
  # Process eggNOG annotations if requested
  if ("eggNOG" %in% include && "eggnog_categories" %in% tables) {
    if (verbose) message("  - Processing eggNOG category annotations...")
    result_data <- .process_eggnog_annotations(con, result_data, base_where, params, verbose)
  }
  
  # Convert to regular data frame and clean up
  result_data <- as.data.frame(result_data)
  result_data$annotation_ids <- NULL  # Remove internal column
  
  if (verbose) {
    message("Processing complete:")
    message("  - Total annotated loci: ", nrow(result_data))
    if ("GO" %in% include) {
      go_count <- sum(!is.na(result_data$go_terms) & result_data$go_terms != "", na.rm = TRUE)
      message("  - Loci with GO annotations: ", go_count)
    }
    if ("KEGG" %in% include) {
      kegg_count <- sum(!is.na(result_data$kegg_pathways) & result_data$kegg_pathways != "", na.rm = TRUE)
      message("  - Loci with KEGG annotations: ", kegg_count)
    }
    if ("Pfam" %in% include) {
      pfam_count <- sum(!is.na(result_data$pfam_domains) & result_data$pfam_domains != "", na.rm = TRUE)
      message("  - Loci with Pfam annotations: ", pfam_count)
    }
    if ("InterPro" %in% include) {
      interpro_count <- sum(!is.na(result_data$interpro_families) & result_data$interpro_families != "", na.rm = TRUE)
      message("  - Loci with InterPro annotations: ", interpro_count)
    }
    if ("eggNOG" %in% include) {
      eggnog_count <- sum(!is.na(result_data$eggnog_categories) & result_data$eggnog_categories != "", na.rm = TRUE)
      message("  - Loci with eggNOG annotations: ", eggnog_count)
    }
  }
  
  # Export CSV if requested
  if (!is.null(export_csv)) {
    if (verbose) message("  - Exporting to CSV: ", export_csv)
    write.csv(result_data, export_csv, row.names = FALSE)
  }
  
  return(result_data)
}

#' Process GO annotations for loci
#' @keywords internal
.process_go_annotations <- function(con, result_data, base_where, params, verbose) {
  
  # Get GO terms for each annotation
  go_query <- paste0("
    SELECT 
      a.annotation_id,
      gt.go_id,
      gt.go_term,
      gt.go_category
    FROM annotations a
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
    JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
    JOIN go_terms gt ON a.annotation_id = gt.annotation_id
    ", base_where
  )
  
  go_data <- if (length(params) > 0) {
    DBI::dbGetQuery(con, go_query, params)
  } else {
    DBI::dbGetQuery(con, go_query)
  }
  
  if (nrow(go_data) > 0) {
    # Aggregate GO terms by annotation
    go_summary <- go_data %>%
      group_by(annotation_id) %>%
      summarise(
        go_terms = paste(unique(go_id), collapse = ";"),
        go_names = paste(unique(go_term), collapse = ";"),
        go_categories = paste(unique(go_category), collapse = ";"),
        .groups = "drop"
      )
    
    # Join with result data through annotation IDs
    result_data <- result_data %>%
      rowwise() %>%
      mutate(
        go_terms = {
          matching_go <- go_summary[go_summary$annotation_id %in% annotation_ids, ]
          if (nrow(matching_go) > 0) {
            paste(unique(unlist(strsplit(matching_go$go_terms, ";"))), collapse = ";")
          } else {
            NA_character_
          }
        },
        go_names = {
          matching_go <- go_summary[go_summary$annotation_id %in% annotation_ids, ]
          if (nrow(matching_go) > 0) {
            paste(unique(unlist(strsplit(matching_go$go_names, ";"))), collapse = ";")
          } else {
            NA_character_
          }
        },
        go_categories = {
          matching_go <- go_summary[go_summary$annotation_id %in% annotation_ids, ]
          if (nrow(matching_go) > 0) {
            paste(unique(unlist(strsplit(matching_go$go_categories, ";"))), collapse = ";")
          } else {
            NA_character_
          }
        }
      ) %>%
      ungroup()
    
    if (verbose) {
      go_count <- sum(!is.na(result_data$go_terms) & result_data$go_terms != "", na.rm = TRUE)
      message("    - ", go_count, " loci have GO annotations")
    }
  } else {
    # Add empty GO columns
    result_data$go_terms <- NA_character_
    result_data$go_names <- NA_character_
    result_data$go_categories <- NA_character_
    if (verbose) message("    - No GO annotations found")
  }
  
  return(result_data)
}

#' Process KEGG annotations with enhancements
#' @keywords internal
.process_kegg_annotations <- function(con, result_data, base_where, params, verbose) {
  
  # Get KEGG references for each annotation
  kegg_query <- paste0("
    SELECT 
      a.annotation_id,
      kr.kegg_id,
      kr.pathway_name
    FROM annotations a
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
    JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
    JOIN kegg_references kr ON a.annotation_id = kr.annotation_id
    ", base_where
  )
  
  kegg_data <- if (length(params) > 0) {
    DBI::dbGetQuery(con, kegg_query, params)
  } else {
    DBI::dbGetQuery(con, kegg_query)
  }
  
  if (nrow(kegg_data) > 0) {
    if (verbose) message("    - Enhancing KEGG pathway names with KEGGREST...")
    
    # Update pathway names using KEGGREST
    kegg_data <- .enhance_kegg_pathway_names(kegg_data, verbose)
    
    if (verbose) message("    - Adding KEGG BRITE classifications...")
    
    # Add KEGG BRITE classifications
    kegg_data <- .add_kegg_brite_classifications(kegg_data, verbose)
    
    if (verbose) message("    - Assigning functional modules...")
    
    # Add functional module assignments
    kegg_data <- .add_kegg_functional_modules(kegg_data, verbose)
    
    # Aggregate KEGG information by annotation
    kegg_summary <- kegg_data %>%
      group_by(annotation_id) %>%
      summarise(
        kegg_pathways = paste(unique(kegg_id), collapse = ";"),
        kegg_pathway_names = paste(unique(pathway_name), collapse = ";"),
        kegg_brite_categories = paste(unique(na.omit(brite_category)), collapse = ";"),
        kegg_modules = paste(unique(na.omit(functional_module)), collapse = ";"),
        .groups = "drop"
      )
    
    # Join with result data through annotation IDs
    result_data <- result_data %>%
      rowwise() %>%
      mutate(
        kegg_pathways = {
          matching_kegg <- kegg_summary[kegg_summary$annotation_id %in% annotation_ids, ]
          if (nrow(matching_kegg) > 0) {
            paste(unique(unlist(strsplit(matching_kegg$kegg_pathways, ";"))), collapse = ";")
          } else {
            NA_character_
          }
        },
        kegg_pathway_names = {
          matching_kegg <- kegg_summary[kegg_summary$annotation_id %in% annotation_ids, ]
          if (nrow(matching_kegg) > 0) {
            paste(unique(unlist(strsplit(matching_kegg$kegg_pathway_names, ";"))), collapse = ";")
          } else {
            NA_character_
          }
        },
        kegg_brite_categories = {
          matching_kegg <- kegg_summary[kegg_summary$annotation_id %in% annotation_ids, ]
          if (nrow(matching_kegg) > 0) {
            cats <- unique(unlist(strsplit(matching_kegg$kegg_brite_categories, ";")))
            cats <- cats[!is.na(cats) & cats != ""]
            if (length(cats) > 0) paste(cats, collapse = ";") else NA_character_
          } else {
            NA_character_
          }
        },
        kegg_modules = {
          matching_kegg <- kegg_summary[kegg_summary$annotation_id %in% annotation_ids, ]
          if (nrow(matching_kegg) > 0) {
            mods <- unique(unlist(strsplit(matching_kegg$kegg_modules, ";")))
            mods <- mods[!is.na(mods) & mods != ""]
            if (length(mods) > 0) paste(mods, collapse = ";") else NA_character_
          } else {
            NA_character_
          }
        }
      ) %>%
      ungroup()
    
    if (verbose) {
      kegg_count <- sum(!is.na(result_data$kegg_pathways) & result_data$kegg_pathways != "", na.rm = TRUE)
      message("    - ", kegg_count, " loci have KEGG annotations")
    }
  } else {
    # Add empty KEGG columns
    result_data$kegg_pathways <- NA_character_
    result_data$kegg_pathway_names <- NA_character_
    result_data$kegg_brite_categories <- NA_character_
    result_data$kegg_modules <- NA_character_
    if (verbose) message("    - No KEGG annotations found")
  }
  
  return(result_data)
}

#' Enhance KEGG pathway names using KEGGREST
#' @keywords internal
.enhance_kegg_pathway_names <- function(kegg_data, verbose) {
  
  # Check if KEGGREST is available
  if (!requireNamespace("KEGGREST", quietly = TRUE)) {
    if (verbose) message("      - KEGGREST not available, using existing pathway names")
    return(kegg_data)
  }
  
  # Find pathways needing name enhancement
  needs_update <- is.na(kegg_data$pathway_name) | 
                  kegg_data$pathway_name == "" | 
                  kegg_data$pathway_name == "-" |
                  grepl("^Gene:", kegg_data$pathway_name)
  
  if (sum(needs_update) == 0) {
    if (verbose) message("      - All pathway names are already descriptive")
    return(kegg_data)
  }
  
  if (verbose) message("      - Updating ", sum(needs_update), " pathway names...")
  
  # Process gene IDs to get pathway information
  gene_ids <- unique(kegg_data$kegg_id[needs_update & grepl("^[a-z]{3}:\\d+", kegg_data$kegg_id)])
  
  for (gene_id in gene_ids[1:min(20, length(gene_ids))]) {  # Limit API calls
    tryCatch({
      gene_info <- KEGGREST::keggGet(gene_id)
      if (length(gene_info) > 0 && !is.null(gene_info[[1]]$PATHWAY)) {
        pathways <- gene_info[[1]]$PATHWAY
        pathway_descriptions <- paste(names(pathways), pathways, sep = ": ")
        
        # Update all matching records
        kegg_data[kegg_data$kegg_id == gene_id & needs_update, "pathway_name"] <- 
          paste(pathway_descriptions, collapse = "; ")
      }
      Sys.sleep(0.2)  # Rate limiting
    }, error = function(e) {
      if (verbose) message("        - Failed to get pathways for ", gene_id)
    })
  }
  
  return(kegg_data)
}

#' Add KEGG BRITE classifications using KEGGREST
#' @keywords internal
.add_kegg_brite_classifications <- function(kegg_data, verbose) {
  
  # For now, use a simplified mapping based on pathway IDs
  # This could be enhanced with actual KEGG BRITE API calls
  
  kegg_data$brite_category <- NA_character_
  
  # Extract pathway IDs from pathway names or kegg_ids
  pathway_patterns <- list(
    "Metabolism" = c("00010", "00020", "00030", "00040", "00051", "00052", "00053", 
                     "00500", "00520", "00620", "00630", "00640", "00650", "00660"),
    "Genetic Information Processing" = c("03010", "03013", "03015", "03018", "03020", 
                                        "03022", "03030", "03040", "03050", "03060"),
    "Environmental Information Processing" = c("02010", "02020", "02024", "02025", 
                                              "02026", "02030", "02040", "02060"),
    "Cellular Processes" = c("04110", "04111", "04112", "04113", "04114", "04120", 
                            "04130", "04136", "04137", "04140", "04141", "04142"),
    "Organismal Systems" = c("04010", "04014", "04015", "04020", "04022", "04024", 
                            "04080", "04142", "04144", "04145", "04146", "04150")
  )
  
  for (category in names(pathway_patterns)) {
    patterns <- pathway_patterns[[category]]
    for (pattern in patterns) {
      matches <- grepl(pattern, kegg_data$kegg_id) | grepl(pattern, kegg_data$pathway_name)
      kegg_data$brite_category[matches & is.na(kegg_data$brite_category)] <- category
    }
  }
  
  return(kegg_data)
}

#' Add functional module assignments
#' @keywords internal  
.add_kegg_functional_modules <- function(kegg_data, verbose) {
  
  # Use simplified functional module classification
  module_patterns <- list(
    "Carbohydrate Metabolism" = c("glycolysis", "gluconeogenesis", "pentose phosphate", 
                                 "citrate cycle", "pyruvate", "starch", "sucrose"),
    "Energy Metabolism" = c("oxidative phosphorylation", "photosynthesis", 
                           "carbon fixation", "nitrogen metabolism"),
    "Lipid Metabolism" = c("fatty acid", "steroid", "glycerolipid", "bile acid"),
    "Amino Acid Metabolism" = c("alanine", "arginine", "aspartate", "glycine", 
                               "leucine", "lysine", "phenylalanine"),
    "Signal Transduction" = c("MAPK", "calcium", "mTOR", "Wnt", "Notch", "insulin"),
    "Cell Cycle" = c("cell cycle", "DNA replication", "p53"),
    "Immune System" = c("complement", "toll-like", "NOD-like", "chemokine")
  )
  
  kegg_data$functional_module <- NA_character_
  
  for (module in names(module_patterns)) {
    patterns <- module_patterns[[module]]
    for (pattern in patterns) {
      matches <- grepl(pattern, kegg_data$pathway_name, ignore.case = TRUE)
      kegg_data$functional_module[matches & is.na(kegg_data$functional_module)] <- module
    }
  }
  
  return(kegg_data)
}

#' Process Pfam annotations for loci
#' @keywords internal
.process_pfam_annotations <- function(con, result_data, base_where, params, verbose) {
  
  # Get Pfam domains for each annotation
  pfam_query <- paste0("
    SELECT 
      a.annotation_id,
      pd.pfam_id,
      pd.domain_name,
      pd.match_status
    FROM annotations a
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
    JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
    JOIN pfam_domains pd ON a.annotation_id = pd.annotation_id
    ", base_where
  )
  
  pfam_data <- if (length(params) > 0) {
    DBI::dbGetQuery(con, pfam_query, params)
  } else {
    DBI::dbGetQuery(con, pfam_query)
  }
  
  if (nrow(pfam_data) > 0) {
    # Aggregate Pfam domains by annotation
    pfam_summary <- pfam_data %>%
      group_by(annotation_id) %>%
      summarise(
        pfam_domains = paste(unique(pfam_id), collapse = ";"),
        pfam_domain_names = paste(unique(na.omit(domain_name)), collapse = ";"),
        .groups = "drop"
      )
    
    # Join with result data through annotation IDs
    result_data <- result_data %>%
      rowwise() %>%
      mutate(
        pfam_domains = {
          matching_pfam <- pfam_summary[pfam_summary$annotation_id %in% annotation_ids, ]
          if (nrow(matching_pfam) > 0) {
            paste(unique(unlist(strsplit(matching_pfam$pfam_domains, ";"))), collapse = ";")
          } else {
            NA_character_
          }
        },
        pfam_domain_names = {
          matching_pfam <- pfam_summary[pfam_summary$annotation_id %in% annotation_ids, ]
          if (nrow(matching_pfam) > 0) {
            names <- unique(unlist(strsplit(matching_pfam$pfam_domain_names, ";")))
            names <- names[!is.na(names) & names != ""]
            if (length(names) > 0) paste(names, collapse = ";") else NA_character_
          } else {
            NA_character_
          }
        }
      ) %>%
      ungroup()
    
    if (verbose) {
      pfam_count <- sum(!is.na(result_data$pfam_domains) & result_data$pfam_domains != "", na.rm = TRUE)
      message("    - ", pfam_count, " loci have Pfam annotations")
    }
  } else {
    # Add empty Pfam columns
    result_data$pfam_domains <- NA_character_
    result_data$pfam_domain_names <- NA_character_
    if (verbose) message("    - No Pfam annotations found")
  }
  
  return(result_data)
}

#' Process InterPro annotations for loci
#' @keywords internal
.process_interpro_annotations <- function(con, result_data, base_where, params, verbose) {
  
  # Get InterPro families for each annotation
  interpro_query <- paste0("
    SELECT 
      a.annotation_id,
      if.interpro_id,
      if.family_name
    FROM annotations a
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
    JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
    JOIN interpro_families if ON a.annotation_id = if.annotation_id
    ", base_where
  )
  
  interpro_data <- if (length(params) > 0) {
    DBI::dbGetQuery(con, interpro_query, params)
  } else {
    DBI::dbGetQuery(con, interpro_query)
  }
  
  if (nrow(interpro_data) > 0) {
    # Aggregate InterPro families by annotation
    interpro_summary <- interpro_data %>%
      group_by(annotation_id) %>%
      summarise(
        interpro_families = paste(unique(interpro_id), collapse = ";"),
        interpro_family_names = paste(unique(na.omit(family_name)), collapse = ";"),
        .groups = "drop"
      )
    
    # Join with result data through annotation IDs
    result_data <- result_data %>%
      rowwise() %>%
      mutate(
        interpro_families = {
          matching_interpro <- interpro_summary[interpro_summary$annotation_id %in% annotation_ids, ]
          if (nrow(matching_interpro) > 0) {
            paste(unique(unlist(strsplit(matching_interpro$interpro_families, ";"))), collapse = ";")
          } else {
            NA_character_
          }
        },
        interpro_family_names = {
          matching_interpro <- interpro_summary[interpro_summary$annotation_id %in% annotation_ids, ]
          if (nrow(matching_interpro) > 0) {
            names <- unique(unlist(strsplit(matching_interpro$interpro_family_names, ";")))
            names <- names[!is.na(names) & names != ""]
            if (length(names) > 0) paste(names, collapse = ";") else NA_character_
          } else {
            NA_character_
          }
        }
      ) %>%
      ungroup()
    
    if (verbose) {
      interpro_count <- sum(!is.na(result_data$interpro_families) & result_data$interpro_families != "", na.rm = TRUE)
      message("    - ", interpro_count, " loci have InterPro annotations")
    }
  } else {
    # Add empty InterPro columns
    result_data$interpro_families <- NA_character_
    result_data$interpro_family_names <- NA_character_
    if (verbose) message("    - No InterPro annotations found")
  }
  
  return(result_data)
}

#' Process eggNOG annotations for loci
#' @keywords internal
.process_eggnog_annotations <- function(con, result_data, base_where, params, verbose) {
  
  # Get eggNOG categories for each annotation
  eggnog_query <- paste0("
    SELECT 
      a.annotation_id,
      ec.eggnog_id,
      ec.taxonomic_scope
    FROM annotations a
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
    JOIN blast_parameters bp ON br.blast_param_id = bp.blast_param_id
    JOIN eggnog_categories ec ON a.annotation_id = ec.annotation_id
    ", base_where
  )
  
  eggnog_data <- if (length(params) > 0) {
    DBI::dbGetQuery(con, eggnog_query, params)
  } else {
    DBI::dbGetQuery(con, eggnog_query)
  }
  
  if (nrow(eggnog_data) > 0) {
    # Aggregate eggNOG categories by annotation
    eggnog_summary <- eggnog_data %>%
      group_by(annotation_id) %>%
      summarise(
        eggnog_categories = paste(unique(eggnog_id), collapse = ";"),
        eggnog_taxonomic_scopes = paste(unique(na.omit(taxonomic_scope)), collapse = ";"),
        .groups = "drop"
      )
    
    # Join with result data through annotation IDs
    result_data <- result_data %>%
      rowwise() %>%
      mutate(
        eggnog_categories = {
          matching_eggnog <- eggnog_summary[eggnog_summary$annotation_id %in% annotation_ids, ]
          if (nrow(matching_eggnog) > 0) {
            paste(unique(unlist(strsplit(matching_eggnog$eggnog_categories, ";"))), collapse = ";")
          } else {
            NA_character_
          }
        },
        eggnog_taxonomic_scopes = {
          matching_eggnog <- eggnog_summary[eggnog_summary$annotation_id %in% annotation_ids, ]
          if (nrow(matching_eggnog) > 0) {
            scopes <- unique(unlist(strsplit(matching_eggnog$eggnog_taxonomic_scopes, ";")))
            scopes <- scopes[!is.na(scopes) & scopes != ""]
            if (length(scopes) > 0) paste(scopes, collapse = ";") else NA_character_
          } else {
            NA_character_
          }
        }
      ) %>%
      ungroup()
    
    if (verbose) {
      eggnog_count <- sum(!is.na(result_data$eggnog_categories) & result_data$eggnog_categories != "", na.rm = TRUE)
      message("    - ", eggnog_count, " loci have eggNOG annotations")
    }
  } else {
    # Add empty eggNOG columns
    result_data$eggnog_categories <- NA_character_
    result_data$eggnog_taxonomic_scopes <- NA_character_
    if (verbose) message("    - No eggNOG annotations found")
  }
  
  return(result_data)
}