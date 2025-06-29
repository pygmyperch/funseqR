#' Descriptive Functional Analysis Functions
#'
#' This file contains functions for descriptive analysis of functional annotations,
#' providing overview and landscape visualization of annotation patterns without
#' statistical testing.
#'

#' Create GO Slim mapping for simplified functional overview
#'
#' Maps detailed GO terms to broader GO Slim categories for high-level functional
#' overview. GO Slim is a subset of GO terms that provides a broad overview of 
#' biological functions without the complexity of the full ontology.
#'
#' @param con Database connection object
#' @param blast_param_id Integer. Optional BLAST parameter ID to filter annotations
#' @param go_slim_set Character. GO Slim set to use: "generic", "plant", "agr", or "pir". Default is "generic"
#' @param include_counts Logical. Include raw counts in output. Default is TRUE
#' @param include_percentages Logical. Include percentages in output. Default is TRUE
#' @param min_count Integer. Minimum count threshold for inclusion. Default is 1
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return Data frame with GO Slim mapping results containing:
#' \itemize{
#'   \item go_slim_category: Broad functional category
#'   \item go_slim_id: GO Slim term ID
#'   \item go_slim_name: GO Slim term name
#'   \item ontology: GO ontology (BP, MF, CC)
#'   \item gene_count: Number of genes annotated (if include_counts = TRUE)
#'   \item percentage: Percentage of total genes (if include_percentages = TRUE)
#'   \item detailed_terms: Number of detailed GO terms mapped to this slim term
#' }
#'
#' @details
#' GO Slim provides a high-level view of functional annotation by grouping
#' detailed GO terms into broader categories. This is particularly useful for:
#' 
#' \itemize{
#'   \item Initial functional overview of datasets
#'   \item Comparing functional profiles between datasets
#'   \item Identifying major functional categories
#'   \item Creating publication-ready functional summaries
#' }
#' 
#' The function maps existing GO terms in the database to their corresponding
#' GO Slim categories using predefined mappings for different GO Slim sets.
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#' 
#' # Basic GO Slim mapping
#' go_slim_summary <- create_go_slim_mapping(con)
#' 
#' # Filtered by BLAST run with specific GO Slim set
#' go_slim_summary <- create_go_slim_mapping(con, blast_param_id = 1, go_slim_set = "agr")
#' 
#' # View results
#' print(go_slim_summary)
#' }
#'
#' @export
create_go_slim_mapping <- function(con, blast_param_id = NULL, go_slim_set = "generic",
                                  include_counts = TRUE, include_percentages = TRUE, 
                                  min_count = 1, verbose = TRUE) {
  
  if (verbose) message("Creating GO Slim mapping...")
  
  # Validate GO Slim set
  valid_sets <- c("generic", "plant", "agr", "pir")
  go_slim_set <- match.arg(go_slim_set, valid_sets)
  
  if (verbose) message("  - Using GO Slim set: ", go_slim_set)
  
  # Check required tables
  tables <- DBI::dbListTables(con)
  required_tables <- c("go_terms", "annotations", "blast_results")
  missing_tables <- required_tables[!required_tables %in% tables]
  
  if (length(missing_tables) > 0) {
    stop("Missing required tables: ", paste(missing_tables, collapse = ", "))
  }
  
  # Build query conditions
  where_conditions <- c()
  params <- list()
  
  if (!is.null(blast_param_id)) {
    where_conditions <- c(where_conditions, "br.blast_param_id = ?")
    params <- c(params, list(blast_param_id))
    if (verbose) message("  - Filtering by BLAST parameter ID: ", blast_param_id)
  }
  
  where_clause <- if (length(where_conditions) > 0) {
    paste("WHERE", paste(where_conditions, collapse = " AND "))
  } else {
    ""
  }
  
  # Get GO terms with annotations
  go_query <- paste0("
    SELECT DISTINCT
      gt.go_id,
      gt.go_term,
      gt.go_category,
      COUNT(DISTINCT a.annotation_id) as annotation_count
    FROM go_terms gt
    JOIN annotations a ON gt.annotation_id = a.annotation_id
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
    ", where_clause, "
    GROUP BY gt.go_id, gt.go_term, gt.go_category
    HAVING COUNT(DISTINCT a.annotation_id) >= ?
    ORDER BY annotation_count DESC
  ")
  
  params <- c(params, list(min_count))
  go_data <- DBI::dbGetQuery(con, go_query, params)
  
  if (nrow(go_data) == 0) {
    if (verbose) message("No GO terms found matching criteria")
    return(data.frame())
  }
  
  if (verbose) message("  - Found ", nrow(go_data), " GO terms to map")
  
  # Apply GO Slim mapping
  go_slim_mapping <- .get_go_slim_mapping(go_slim_set)
  
  # Map GO terms to GO Slim categories
  mapped_data <- .map_go_terms_to_slim(go_data, go_slim_mapping, verbose)
  
  # Calculate summary statistics
  slim_summary <- mapped_data %>%
    dplyr::group_by(go_slim_id, go_slim_name, go_slim_category, ontology) %>%
    dplyr::summarise(
      gene_count = sum(annotation_count),
      detailed_terms = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::arrange(desc(gene_count))
  
  # Add percentages if requested
  if (include_percentages) {
    total_genes <- sum(slim_summary$gene_count)
    slim_summary$percentage <- round((slim_summary$gene_count / total_genes) * 100, 2)
  }
  
  # Filter by minimum counts if requested
  if (!include_counts) {
    slim_summary$gene_count <- NULL
  }
  
  if (verbose) {
    message("GO Slim mapping complete:")
    message("  - Mapped to ", nrow(slim_summary), " GO Slim categories")
    message("  - Total genes: ", sum(slim_summary$gene_count))
  }
  
  return(as.data.frame(slim_summary))
}

#' Analyze KEGG BRITE hierarchy for pathway organization
#'
#' Provides hierarchical analysis of KEGG pathways using the BRITE functional
#' classification system. BRITE organizes KEGG pathways into hierarchical
#' categories for systematic functional analysis.
#'
#' @param con Database connection object
#' @param blast_param_id Integer. Optional BLAST parameter ID to filter annotations
#' @param hierarchy_level Integer. BRITE hierarchy level to analyze (1-3). Default is 2
#' @param include_pathway_details Logical. Include individual pathway information. Default is TRUE
#' @param min_count Integer. Minimum pathway count for inclusion. Default is 1
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return List containing:
#' \itemize{
#'   \item summary: Data frame with BRITE category summaries
#'   \item hierarchy: Complete hierarchical structure
#'   \item pathways: Individual pathway details (if include_pathway_details = TRUE)
#' }
#'
#' @details
#' KEGG BRITE provides a hierarchical classification of biological systems:
#' 
#' \strong{Level 1}: Major functional categories
#' \itemize{
#'   \item Metabolism
#'   \item Genetic Information Processing  
#'   \item Environmental Information Processing
#'   \item Cellular Processes
#'   \item Organismal Systems
#' }
#' 
#' \strong{Level 2}: Functional subcategories within each major category
#' 
#' \strong{Level 3}: Specific pathway groups
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#' 
#' # Analyze BRITE hierarchy at level 2
#' brite_analysis <- analyze_kegg_brite_hierarchy(con, hierarchy_level = 2)
#' 
#' # View summary
#' print(brite_analysis$summary)
#' 
#' # View hierarchy structure
#' print(brite_analysis$hierarchy)
#' }
#'
#' @export
analyze_kegg_brite_hierarchy <- function(con, blast_param_id = NULL, hierarchy_level = 2,
                                        include_pathway_details = TRUE, min_count = 1, verbose = TRUE) {
  
  if (verbose) message("Analyzing KEGG BRITE hierarchy...")
  
  # Validate hierarchy level
  if (!hierarchy_level %in% 1:3) {
    stop("hierarchy_level must be 1, 2, or 3")
  }
  
  if (verbose) message("  - Analyzing hierarchy level: ", hierarchy_level)
  
  # Check required tables
  tables <- DBI::dbListTables(con)
  required_tables <- c("kegg_references", "annotations", "blast_results")
  missing_tables <- required_tables[!required_tables %in% tables]
  
  if (length(missing_tables) > 0) {
    stop("Missing required tables: ", paste(missing_tables, collapse = ", "))
  }
  
  # Build query conditions
  where_conditions <- c()
  params <- list()
  
  if (!is.null(blast_param_id)) {
    where_conditions <- c(where_conditions, "br.blast_param_id = ?")
    params <- c(params, list(blast_param_id))
    if (verbose) message("  - Filtering by BLAST parameter ID: ", blast_param_id)
  }
  
  where_clause <- if (length(where_conditions) > 0) {
    paste("WHERE", paste(where_conditions, collapse = " AND "))
  } else {
    ""
  }
  
  # Get KEGG pathways with annotations
  kegg_query <- paste0("
    SELECT DISTINCT
      kr.kegg_id,
      kr.pathway_name,
      COUNT(DISTINCT a.annotation_id) as annotation_count
    FROM kegg_references kr
    JOIN annotations a ON kr.annotation_id = a.annotation_id
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
    ", where_clause, "
    GROUP BY kr.kegg_id, kr.pathway_name
    HAVING COUNT(DISTINCT a.annotation_id) >= ?
    ORDER BY annotation_count DESC
  ")
  
  params <- c(params, list(min_count))
  kegg_data <- DBI::dbGetQuery(con, kegg_query, params)
  
  if (nrow(kegg_data) == 0) {
    if (verbose) message("No KEGG pathways found matching criteria")
    return(list(summary = data.frame(), hierarchy = list(), pathways = data.frame()))
  }
  
  if (verbose) message("  - Found ", nrow(kegg_data), " KEGG pathways to analyze")
  
  # Get BRITE hierarchy mapping
  brite_hierarchy <- .get_kegg_brite_hierarchy()
  
  # Map pathways to BRITE hierarchy
  mapped_pathways <- .map_kegg_to_brite(kegg_data, brite_hierarchy, hierarchy_level, verbose)
  
  # Create summary at requested hierarchy level
  hierarchy_summary <- mapped_pathways %>%
    dplyr::group_by(brite_level1, brite_level2, brite_level3) %>%
    dplyr::summarise(
      pathway_count = dplyr::n(),
      total_genes = sum(annotation_count),
      .groups = "drop"
    ) %>%
    dplyr::arrange(desc(total_genes))
  
  # Filter by hierarchy level
  if (hierarchy_level == 1) {
    summary_data <- hierarchy_summary %>%
      dplyr::group_by(brite_category = brite_level1) %>%
      dplyr::summarise(
        pathway_count = sum(pathway_count),
        total_genes = sum(total_genes),
        .groups = "drop"
      )
  } else if (hierarchy_level == 2) {
    summary_data <- hierarchy_summary %>%
      dplyr::group_by(brite_level1, brite_category = brite_level2) %>%
      dplyr::summarise(
        pathway_count = sum(pathway_count),
        total_genes = sum(total_genes),
        .groups = "drop"
      )
  } else {
    summary_data <- hierarchy_summary %>%
      dplyr::mutate(brite_category = brite_level3) %>%
      dplyr::select(brite_level1, brite_level2, brite_category, pathway_count, total_genes)
  }
  
  # Add percentages
  summary_data$percentage <- round((summary_data$total_genes / sum(summary_data$total_genes)) * 100, 2)
  
  # Prepare result
  result <- list(
    summary = as.data.frame(summary_data),
    hierarchy = as.list(hierarchy_summary)
  )
  
  if (include_pathway_details) {
    result$pathways <- as.data.frame(mapped_pathways)
  }
  
  if (verbose) {
    message("KEGG BRITE analysis complete:")
    message("  - Analyzed ", nrow(summary_data), " BRITE categories")
    message("  - Total pathways: ", sum(summary_data$pathway_count))
    message("  - Total genes: ", sum(summary_data$total_genes))
  }
  
  return(result)
}

#' Create functional landscape visualization
#'
#' Generates comprehensive functional landscape visualization combining multiple
#' annotation types (GO, KEGG, Pfam, InterPro, eggNOG) to provide an integrated
#' view of functional annotation patterns.
#'
#' @param con Database connection object
#' @param blast_param_id Integer. Optional BLAST parameter ID to filter annotations
#' @param annotation_types Character vector. Annotation types to include. Default is c("GO", "KEGG", "Pfam", "InterPro", "eggNOG")
#' @param plot_type Character. Type of visualization: "treemap", "sunburst", "network", or "heatmap". Default is "treemap"
#' @param top_n Integer. Number of top terms per annotation type to include. Default is 10
#' @param min_genes Integer. Minimum gene count for inclusion. Default is 2
#' @param color_palette Character. Color palette to use. Default is "Set3"
#' @param export_data Logical. Export underlying data. Default is FALSE
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return List containing:
#' \itemize{
#'   \item plot: ggplot object with functional landscape visualization
#'   \item summary: Summary statistics across annotation types
#'   \item data: Raw data used for visualization (if export_data = TRUE)
#' }
#'
#' @details
#' The functional landscape provides an integrated view of annotations by:
#' 
#' \itemize{
#'   \item Combining multiple annotation databases
#'   \item Showing relative importance of functional categories
#'   \item Highlighting annotation coverage and patterns
#'   \item Enabling comparison between datasets
#' }
#' 
#' Different plot types offer different perspectives:
#' \itemize{
#'   \item \strong{treemap}: Hierarchical area-based visualization
#'   \item \strong{sunburst}: Circular hierarchical visualization
#'   \item \strong{network}: Network-based functional relationships
#'   \item \strong{heatmap}: Matrix-based comparison view
#' }
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#' 
#' # Create treemap functional landscape
#' landscape <- create_functional_landscape(con, plot_type = "treemap")
#' print(landscape$plot)
#' 
#' # Create landscape with specific annotation types
#' landscape <- create_functional_landscape(con, 
#'                                         annotation_types = c("GO", "Pfam", "InterPro"),
#'                                         plot_type = "sunburst",
#'                                         top_n = 15)
#' print(landscape$summary)
#' }
#'
#' @export
create_functional_landscape <- function(con, blast_param_id = NULL, 
                                       annotation_types = c("GO", "KEGG", "Pfam", "InterPro", "eggNOG"),
                                       plot_type = "treemap", top_n = 10, min_genes = 2,
                                       color_palette = "Set3", export_data = FALSE, verbose = TRUE) {
  
  if (verbose) message("Creating functional landscape visualization...")
  
  # Validate inputs
  valid_types <- c("GO", "KEGG", "Pfam", "InterPro", "eggNOG")
  annotation_types <- match.arg(annotation_types, valid_types, several.ok = TRUE)
  
  valid_plots <- c("treemap", "sunburst", "network", "heatmap")
  plot_type <- match.arg(plot_type, valid_plots)
  
  if (verbose) {
    message("  - Annotation types: ", paste(annotation_types, collapse = ", "))
    message("  - Plot type: ", plot_type)
    message("  - Top terms per type: ", top_n)
  }
  
  # Collect data from each annotation type
  landscape_data <- list()
  
  for (annotation_type in annotation_types) {
    if (verbose) message("  - Processing ", annotation_type, " annotations...")
    
    type_data <- .get_annotation_landscape_data(con, annotation_type, blast_param_id, 
                                               top_n, min_genes, verbose)
    
    if (nrow(type_data) > 0) {
      landscape_data[[annotation_type]] <- type_data
    }
  }
  
  if (length(landscape_data) == 0) {
    stop("No annotation data found for any requested types")
  }
  
  # Combine data from all annotation types
  combined_data <- do.call(rbind, lapply(names(landscape_data), function(type) {
    data <- landscape_data[[type]]
    data$annotation_type <- type
    return(data)
  }))
  
  # Create visualization based on plot type
  plot_obj <- switch(plot_type,
    "treemap" = .create_treemap_landscape(combined_data, color_palette),
    "sunburst" = .create_sunburst_landscape(combined_data, color_palette),
    "network" = .create_network_landscape(combined_data, color_palette),
    "heatmap" = .create_heatmap_landscape(combined_data, color_palette)
  )
  
  # Create summary statistics
  summary_stats <- combined_data %>%
    dplyr::group_by(annotation_type) %>%
    dplyr::summarise(
      unique_terms = dplyr::n(),
      total_genes = sum(gene_count),
      avg_genes_per_term = round(mean(gene_count), 1),
      max_genes_per_term = max(gene_count),
      .groups = "drop"
    ) %>%
    dplyr::arrange(desc(total_genes))
  
  # Calculate coverage percentages
  total_unique_genes <- length(unique(unlist(strsplit(combined_data$gene_list, ";"))))
  summary_stats$coverage_percentage <- round((summary_stats$total_genes / total_unique_genes) * 100, 1)
  
  # Prepare result
  result <- list(
    plot = plot_obj,
    summary = as.data.frame(summary_stats)
  )
  
  if (export_data) {
    result$data <- as.data.frame(combined_data)
  }
  
  if (verbose) {
    message("Functional landscape complete:")
    message("  - Total annotation types: ", length(landscape_data))
    message("  - Total functional terms: ", nrow(combined_data))
    message("  - Total unique genes: ", total_unique_genes)
  }
  
  return(result)
}

# Helper functions for descriptive analysis

#' Get GO Slim mapping data
#' @keywords internal
.get_go_slim_mapping <- function(go_slim_set) {
  
  # Define GO Slim mappings for different sets
  # This is a simplified mapping - in practice, these would be loaded from official GO Slim files
  
  if (go_slim_set == "generic") {
    # Generic GO Slim - broad functional categories
    mapping <- list(
      # Biological Process
      "GO:0008150" = list(name = "biological_process", category = "root", ontology = "BP"),
      "GO:0009987" = list(name = "cellular process", category = "cellular", ontology = "BP"),
      "GO:0008152" = list(name = "metabolic process", category = "metabolism", ontology = "BP"),
      "GO:0050896" = list(name = "response to stimulus", category = "response", ontology = "BP"),
      "GO:0065007" = list(name = "biological regulation", category = "regulation", ontology = "BP"),
      "GO:0051179" = list(name = "localization", category = "transport", ontology = "BP"),
      "GO:0032501" = list(name = "multicellular organismal process", category = "development", ontology = "BP"),
      "GO:0032502" = list(name = "developmental process", category = "development", ontology = "BP"),
      "GO:0040011" = list(name = "locomotion", category = "behavior", ontology = "BP"),
      "GO:0071840" = list(name = "cellular component organization", category = "organization", ontology = "BP"),
      
      # Molecular Function
      "GO:0003674" = list(name = "molecular_function", category = "root", ontology = "MF"),
      "GO:0003824" = list(name = "catalytic activity", category = "catalytic", ontology = "MF"),
      "GO:0005488" = list(name = "binding", category = "binding", ontology = "MF"),
      "GO:0005215" = list(name = "transporter activity", category = "transport", ontology = "MF"),
      "GO:0038023" = list(name = "signaling receptor activity", category = "signaling", ontology = "MF"),
      "GO:0005198" = list(name = "structural molecule activity", category = "structural", ontology = "MF"),
      "GO:0030234" = list(name = "enzyme regulator activity", category = "regulation", ontology = "MF"),
      "GO:0140110" = list(name = "transcription regulator activity", category = "regulation", ontology = "MF"),
      
      # Cellular Component
      "GO:0005575" = list(name = "cellular_component", category = "root", ontology = "CC"),
      "GO:0005623" = list(name = "cell", category = "cell", ontology = "CC"),
      "GO:0032991" = list(name = "protein-containing complex", category = "complex", ontology = "CC"),
      "GO:0005576" = list(name = "extracellular region", category = "extracellular", ontology = "CC"),
      "GO:0016020" = list(name = "membrane", category = "membrane", ontology = "CC"),
      "GO:0030054" = list(name = "cell junction", category = "junction", ontology = "CC")
    )
  } else {
    # For other GO Slim sets, use a simplified version
    # In practice, these would be loaded from appropriate GO Slim files
    mapping <- list(
      # Basic categories for any GO Slim set
      "GO:0008150" = list(name = "biological_process", category = "process", ontology = "BP"),
      "GO:0003674" = list(name = "molecular_function", category = "function", ontology = "MF"),
      "GO:0005575" = list(name = "cellular_component", category = "component", ontology = "CC")
    )
  }
  
  return(mapping)
}

#' Map GO terms to GO Slim categories
#' @keywords internal
.map_go_terms_to_slim <- function(go_data, go_slim_mapping, verbose) {
  
  if (verbose) message("    - Mapping GO terms to slim categories...")
  
  # Initialize result data frame
  mapped_data <- data.frame(
    go_id = character(),
    go_term = character(),
    go_category = character(),
    annotation_count = integer(),
    go_slim_id = character(),
    go_slim_name = character(),
    go_slim_category = character(),
    ontology = character(),
    stringsAsFactors = FALSE
  )
  
  # Simple mapping approach - in practice this would use GO hierarchy
  for (i in 1:nrow(go_data)) {
    go_id <- go_data$go_id[i]
    go_term <- go_data$go_term[i]
    go_category <- go_data$go_category[i]
    count <- go_data$annotation_count[i]
    
    # Find best matching GO Slim category based on ontology and keywords
    slim_match <- .find_best_slim_match(go_term, go_category, go_slim_mapping)
    
    if (!is.null(slim_match)) {
      mapped_data <- rbind(mapped_data, data.frame(
        go_id = go_id,
        go_term = go_term,
        go_category = go_category,
        annotation_count = count,
        go_slim_id = slim_match$id,
        go_slim_name = slim_match$name,
        go_slim_category = slim_match$category,
        ontology = slim_match$ontology,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  if (verbose) message("    - Mapped ", nrow(mapped_data), " terms to slim categories")
  
  return(mapped_data)
}

#' Find best GO Slim match for a GO term
#' @keywords internal
.find_best_slim_match <- function(go_term, go_category, go_slim_mapping) {
  
  # Convert GO category to ontology code
  ontology <- switch(go_category,
    "biological_process" = "BP",
    "molecular_function" = "MF", 
    "cellular_component" = "CC",
    go_category  # fallback
  )
  
  # Filter mapping by ontology
  ontology_mapping <- go_slim_mapping[sapply(go_slim_mapping, function(x) x$ontology == ontology)]
  
  if (length(ontology_mapping) == 0) {
    return(NULL)
  }
  
  # Simple keyword-based matching
  go_term_lower <- tolower(go_term)
  
  # Define keyword patterns for matching
  keyword_patterns <- list(
    "metabolic" = c("metabolic", "metabolism", "biosynthetic", "catabolic"),
    "cellular" = c("cellular", "cell"),
    "transport" = c("transport", "localization", "secretion"),
    "response" = c("response", "stimulus", "stress"),
    "regulation" = c("regulation", "control", "modulation"),
    "development" = c("development", "differentiation", "morphogenesis"),
    "catalytic" = c("catalytic", "enzyme", "activity"),
    "binding" = c("binding", "bind"),
    "membrane" = c("membrane", "transmembrane"),
    "nucleus" = c("nucleus", "nuclear"),
    "extracellular" = c("extracellular", "secreted")
  )
  
  # Try to match based on keywords
  for (slim_id in names(ontology_mapping)) {
    slim_info <- ontology_mapping[[slim_id]]
    slim_category <- slim_info$category
    
    if (slim_category %in% names(keyword_patterns)) {
      patterns <- keyword_patterns[[slim_category]]
      if (any(sapply(patterns, function(p) grepl(p, go_term_lower)))) {
        return(list(
          id = slim_id,
          name = slim_info$name,
          category = slim_info$category,
          ontology = slim_info$ontology
        ))
      }
    }
  }
  
  # Default to root term for the ontology
  root_terms <- list(
    "BP" = "GO:0008150",
    "MF" = "GO:0003674", 
    "CC" = "GO:0005575"
  )
  
  root_id <- root_terms[[ontology]]
  if (root_id %in% names(ontology_mapping)) {
    root_info <- ontology_mapping[[root_id]]
    return(list(
      id = root_id,
      name = root_info$name,
      category = root_info$category,
      ontology = root_info$ontology
    ))
  }
  
  return(NULL)
}

#' Get KEGG BRITE hierarchy data
#' @keywords internal
.get_kegg_brite_hierarchy <- function() {
  
  # KEGG BRITE functional hierarchy
  # This is a simplified version - in practice would be loaded from KEGG BRITE files
  
  hierarchy <- list(
    # Level 1: Major categories
    "Metabolism" = list(
      "Carbohydrate metabolism" = c("00010", "00020", "00030", "00040", "00051", "00052", "00053"),
      "Energy metabolism" = c("00190", "00195", "00196", "00710", "00720"),
      "Lipid metabolism" = c("00061", "00062", "00071", "00072", "00073", "00100"),
      "Nucleotide metabolism" = c("00230", "00240"),
      "Amino acid metabolism" = c("00250", "00260", "00270", "00280", "00290")
    ),
    
    "Genetic Information Processing" = list(
      "Transcription" = c("03020", "03022", "03040"),
      "Translation" = c("03010", "03013", "03015"),
      "Folding and degradation" = c("03060", "04141"),
      "Replication and repair" = c("03030", "03410", "03420", "03430")
    ),
    
    "Environmental Information Processing" = list(
      "Membrane transport" = c("02010", "02020"),
      "Signal transduction" = c("02024", "02025", "02026", "04010", "04014"),
      "Signaling pathways" = c("04080", "04142", "04144", "04145")
    ),
    
    "Cellular Processes" = list(
      "Cell growth and death" = c("04110", "04111", "04112", "04113", "04114", "04115"),
      "Cellular community" = c("04120", "04130", "04136", "04137"),
      "Cell motility" = c("04810", "04814")
    ),
    
    "Organismal Systems" = list(
      "Immune system" = c("04610", "04611", "04612", "04620", "04621", "04622"),
      "Endocrine system" = c("04910", "04911", "04912", "04913", "04914", "04915"),
      "Nervous system" = c("04720", "04721", "04722", "04723", "04724", "04725")
    )
  )
  
  return(hierarchy)
}

#' Map KEGG pathways to BRITE hierarchy
#' @keywords internal
.map_kegg_to_brite <- function(kegg_data, brite_hierarchy, hierarchy_level, verbose) {
  
  if (verbose) message("    - Mapping KEGG pathways to BRITE hierarchy...")
  
  # Initialize result data frame
  mapped_data <- data.frame(
    kegg_id = character(),
    pathway_name = character(),
    annotation_count = integer(),
    brite_level1 = character(),
    brite_level2 = character(),
    brite_level3 = character(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:nrow(kegg_data)) {
    kegg_id <- kegg_data$kegg_id[i]
    pathway_name <- kegg_data$pathway_name[i]
    count <- kegg_data$annotation_count[i]
    
    # Extract pathway number from KEGG ID
    pathway_num <- .extract_pathway_number(kegg_id)
    
    # Find BRITE classification
    brite_class <- .find_brite_classification(pathway_num, brite_hierarchy)
    
    if (!is.null(brite_class)) {
      mapped_data <- rbind(mapped_data, data.frame(
        kegg_id = kegg_id,
        pathway_name = pathway_name,
        annotation_count = count,
        brite_level1 = brite_class$level1,
        brite_level2 = brite_class$level2,
        brite_level3 = ifelse(is.null(brite_class$level3), NA_character_, brite_class$level3),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  if (verbose) message("    - Mapped ", nrow(mapped_data), " pathways to BRITE categories")
  
  return(mapped_data)
}

#' Extract pathway number from KEGG ID
#' @keywords internal
.extract_pathway_number <- function(kegg_id) {
  # Extract pathway number (e.g., "00010" from "map00010" or "hsa00010")
  if (grepl("\\d{5}", kegg_id)) {
    # Use base R instead of stringr for better compatibility
    match <- regmatches(kegg_id, regexpr("\\d{5}", kegg_id))
    return(if(length(match) > 0) match[1] else NULL)
  }
  return(NULL)
}

#' Find BRITE classification for pathway
#' @keywords internal
.find_brite_classification <- function(pathway_num, brite_hierarchy) {
  
  if (is.null(pathway_num)) {
    return(NULL)
  }
  
  for (level1 in names(brite_hierarchy)) {
    level2_categories <- brite_hierarchy[[level1]]
    
    for (level2 in names(level2_categories)) {
      pathway_list <- level2_categories[[level2]]
      
      if (pathway_num %in% pathway_list) {
        return(list(
          level1 = level1,
          level2 = level2,
          level3 = NULL  # Could be extended for level 3
        ))
      }
    }
  }
  
  # Default classification if not found
  return(list(
    level1 = "Unknown",
    level2 = "Unclassified",
    level3 = NULL
  ))
}

#' Get annotation landscape data for visualization
#' @keywords internal
.get_annotation_landscape_data <- function(con, annotation_type, blast_param_id, top_n, min_genes, verbose) {
  
  # Build query based on annotation type
  if (annotation_type == "GO") {
    query <- .build_go_landscape_query(blast_param_id)
  } else if (annotation_type == "KEGG") {
    query <- .build_kegg_landscape_query(blast_param_id)
  } else if (annotation_type == "Pfam") {
    query <- .build_pfam_landscape_query(blast_param_id)
  } else if (annotation_type == "InterPro") {
    query <- .build_interpro_landscape_query(blast_param_id)
  } else if (annotation_type == "eggNOG") {
    query <- .build_eggnog_landscape_query(blast_param_id)
  } else {
    return(data.frame())
  }
  
  # Execute query
  params <- list()
  if (!is.null(blast_param_id)) {
    params <- list(blast_param_id, min_genes)
  } else {
    params <- list(min_genes)
  }
  
  raw_data <- DBI::dbGetQuery(con, query, params)
  
  if (nrow(raw_data) == 0) {
    return(data.frame())
  }
  
  # Process and format data
  processed_data <- raw_data %>%
    dplyr::arrange(desc(gene_count)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(
      gene_list = paste(unique_genes, collapse = ";"),
      percentage = round((gene_count / sum(gene_count)) * 100, 1)
    )
  
  return(processed_data)
}

#' Build GO landscape query
#' @keywords internal
.build_go_landscape_query <- function(blast_param_id) {
  
  base_query <- "
    SELECT 
      gt.go_id as term_id,
      gt.go_term as term_name,
      gt.go_category as category,
      COUNT(DISTINCT a.annotation_id) as gene_count,
      GROUP_CONCAT(DISTINCT a.uniprot_accession) as unique_genes
    FROM go_terms gt
    JOIN annotations a ON gt.annotation_id = a.annotation_id
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
  "
  
  if (!is.null(blast_param_id)) {
    base_query <- paste0(base_query, "
    WHERE br.blast_param_id = ?
    GROUP BY gt.go_id, gt.go_term, gt.go_category
    HAVING COUNT(DISTINCT a.annotation_id) >= ?
    ORDER BY gene_count DESC")
  } else {
    base_query <- paste0(base_query, "
    GROUP BY gt.go_id, gt.go_term, gt.go_category
    HAVING COUNT(DISTINCT a.annotation_id) >= ?
    ORDER BY gene_count DESC")
  }
  
  return(base_query)
}

#' Build KEGG landscape query
#' @keywords internal
.build_kegg_landscape_query <- function(blast_param_id) {
  
  base_query <- "
    SELECT 
      kr.kegg_id as term_id,
      kr.pathway_name as term_name,
      'PATHWAY' as category,
      COUNT(DISTINCT a.annotation_id) as gene_count,
      GROUP_CONCAT(DISTINCT a.uniprot_accession) as unique_genes
    FROM kegg_references kr
    JOIN annotations a ON kr.annotation_id = a.annotation_id
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
  "
  
  if (!is.null(blast_param_id)) {
    base_query <- paste0(base_query, "
    WHERE br.blast_param_id = ?
    GROUP BY kr.kegg_id, kr.pathway_name
    HAVING COUNT(DISTINCT a.annotation_id) >= ?
    ORDER BY gene_count DESC")
  } else {
    base_query <- paste0(base_query, "
    GROUP BY kr.kegg_id, kr.pathway_name
    HAVING COUNT(DISTINCT a.annotation_id) >= ?
    ORDER BY gene_count DESC")
  }
  
  return(base_query)
}

#' Build Pfam landscape query
#' @keywords internal
.build_pfam_landscape_query <- function(blast_param_id) {
  
  base_query <- "
    SELECT 
      pd.pfam_id as term_id,
      COALESCE(pd.domain_name, pd.pfam_id) as term_name,
      'DOMAIN' as category,
      COUNT(DISTINCT a.annotation_id) as gene_count,
      GROUP_CONCAT(DISTINCT a.uniprot_accession) as unique_genes
    FROM pfam_domains pd
    JOIN annotations a ON pd.annotation_id = a.annotation_id
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
  "
  
  if (!is.null(blast_param_id)) {
    base_query <- paste0(base_query, "
    WHERE br.blast_param_id = ?
    GROUP BY pd.pfam_id, pd.domain_name
    HAVING COUNT(DISTINCT a.annotation_id) >= ?
    ORDER BY gene_count DESC")
  } else {
    base_query <- paste0(base_query, "
    GROUP BY pd.pfam_id, pd.domain_name
    HAVING COUNT(DISTINCT a.annotation_id) >= ?
    ORDER BY gene_count DESC")
  }
  
  return(base_query)
}

#' Build InterPro landscape query
#' @keywords internal
.build_interpro_landscape_query <- function(blast_param_id) {
  
  base_query <- "
    SELECT 
      if.interpro_id as term_id,
      COALESCE(if.family_name, if.interpro_id) as term_name,
      'FAMILY' as category,
      COUNT(DISTINCT a.annotation_id) as gene_count,
      GROUP_CONCAT(DISTINCT a.uniprot_accession) as unique_genes
    FROM interpro_families if
    JOIN annotations a ON if.annotation_id = a.annotation_id
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
  "
  
  if (!is.null(blast_param_id)) {
    base_query <- paste0(base_query, "
    WHERE br.blast_param_id = ?
    GROUP BY if.interpro_id, if.family_name
    HAVING COUNT(DISTINCT a.annotation_id) >= ?
    ORDER BY gene_count DESC")
  } else {
    base_query <- paste0(base_query, "
    GROUP BY if.interpro_id, if.family_name
    HAVING COUNT(DISTINCT a.annotation_id) >= ?
    ORDER BY gene_count DESC")
  }
  
  return(base_query)
}

#' Build eggNOG landscape query
#' @keywords internal
.build_eggnog_landscape_query <- function(blast_param_id) {
  
  base_query <- "
    SELECT 
      ec.eggnog_id as term_id,
      COALESCE(ec.taxonomic_scope, ec.eggnog_id) as term_name,
      'ORTHOLOGY' as category,
      COUNT(DISTINCT a.annotation_id) as gene_count,
      GROUP_CONCAT(DISTINCT a.uniprot_accession) as unique_genes
    FROM eggnog_categories ec
    JOIN annotations a ON ec.annotation_id = a.annotation_id
    JOIN blast_results br ON a.blast_result_id = br.blast_result_id
  "
  
  if (!is.null(blast_param_id)) {
    base_query <- paste0(base_query, "
    WHERE br.blast_param_id = ?
    GROUP BY ec.eggnog_id, ec.taxonomic_scope
    HAVING COUNT(DISTINCT a.annotation_id) >= ?
    ORDER BY gene_count DESC")
  } else {
    base_query <- paste0(base_query, "
    GROUP BY ec.eggnog_id, ec.taxonomic_scope
    HAVING COUNT(DISTINCT a.annotation_id) >= ?
    ORDER BY gene_count DESC")
  }
  
  return(base_query)
}

# Visualization helper functions

#' Create treemap landscape visualization
#' @keywords internal
.create_treemap_landscape <- function(combined_data, color_palette) {
  
  if (!requireNamespace("treemapify", quietly = TRUE)) {
    stop("treemapify package is required for treemap visualization")
  }
  
  # Prepare data for treemap
  treemap_data <- combined_data %>%
    dplyr::mutate(
      label = paste0(term_name, "\n(", gene_count, ")"),
      size = gene_count
    )
  
  # Create treemap
  p <- ggplot2::ggplot(treemap_data, ggplot2::aes(area = size, fill = annotation_type, label = label)) +
    treemapify::geom_treemap() +
    treemapify::geom_treemap_text(colour = "white", place = "centre", size = 8) +
    ggplot2::scale_fill_brewer(palette = color_palette, name = "Annotation Type") +
    ggplot2::theme_void() +
    ggplot2::labs(title = "Functional Landscape - Treemap View",
                  subtitle = "Size represents number of genes annotated")
  
  return(p)
}

#' Create sunburst landscape visualization (placeholder)
#' @keywords internal
.create_sunburst_landscape <- function(combined_data, color_palette) {
  
  # Placeholder for sunburst plot - would require additional packages like plotly or sunburstR
  message("Sunburst visualization not yet implemented - showing treemap instead")
  return(.create_treemap_landscape(combined_data, color_palette))
}

#' Create network landscape visualization (placeholder)
#' @keywords internal
.create_network_landscape <- function(combined_data, color_palette) {
  
  # Placeholder for network plot - would require packages like igraph or ggraph
  message("Network visualization not yet implemented - showing heatmap instead")
  return(.create_heatmap_landscape(combined_data, color_palette))
}

#' Create heatmap landscape visualization
#' @keywords internal
.create_heatmap_landscape <- function(combined_data, color_palette) {
  
  # Create heatmap showing annotation types vs top terms
  heatmap_data <- combined_data %>%
    dplyr::group_by(annotation_type) %>%
    dplyr::slice_head(n = 10) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      term_short = ifelse(nchar(term_name) > 40, paste0(substr(term_name, 1, 37), "..."), term_name),
      log_count = log10(gene_count + 1)
    )
  
  p <- ggplot2::ggplot(heatmap_data, ggplot2::aes(x = annotation_type, y = reorder(term_short, gene_count))) +
    ggplot2::geom_tile(ggplot2::aes(fill = log_count), color = "white") +
    ggplot2::scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Log10\n(Gene Count)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 8)) +
    ggplot2::labs(title = "Functional Landscape - Heatmap View",
                  subtitle = "Top functional terms by annotation type",
                  x = "Annotation Type",
                  y = "Functional Terms")
  
  return(p)
}