#' Project Management Functions for funseqR
#'
#' Functions for managing funseqR projects, including data export and project duplication.
#'

#' Export data from existing project to create new project
#'
#' Creates a new funseqR database by copying data from an existing project.
#' This is a simple, efficient bulk copy operation that preserves all IDs and relationships.
#'
#' @param source_con Database connection object for the source project
#' @param target_db_path Character string specifying the path for the new database file
#' @param tables Character vector specifying which tables to export. If NULL, exports all available tables.
#'   Default table order respects foreign key dependencies.
#' @param force Logical. If TRUE, overwrite existing database file. Default is FALSE
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return Database connection object to the newly created project database
#'
#' @details
#' This function performs a simple bulk copy operation to duplicate funseqR project data.
#' Unlike complex export solutions, this preserves all primary keys and foreign key relationships
#' by copying data exactly as stored in the source database.
#'
#' \strong{Default Export Tables (in dependency order):}
#' \itemize{
#'   \item input_files - File registrations
#'   \item reference_genomes - Reference genome metadata  
#'   \item reference_sequences - Reference genome sequences
#'   \item vcf_data - Variant data
#'   \item candidate_loci - Candidate loci coordinates
#'   \item locus_statistics - Statistical values for loci
#'   \item flanking_sequences - Flanking sequences around variants
#'   \item blast_parameters - BLAST search parameters
#'   \item blast_results - BLAST search results
#'   \item annotations - Functional annotations
#'   \item go_terms - Gene Ontology annotations
#'   \item kegg_references - KEGG pathway annotations
#'   \item pfam_domains - Pfam domain annotations
#'   \item interpro_families - InterPro family annotations
#'   \item eggnog_categories - eggNOG category annotations
#'   \item uniprot_cache - UniProt data cache
#' }
#'
#' @examples
#' \dontrun{
#' # Connect to existing project
#' source_con <- connect_funseq_db("original_project.db")
#' 
#' # Complete project copy
#' new_con <- export_project_data(source_con, "complete_copy.db")
#' 
#' # Copy only essential data
#' new_con <- export_project_data(
#'   source_con, "essential.db",
#'   tables = c("input_files", "reference_genomes", "reference_sequences", "vcf_data")
#' )
#' 
#' # Copy up to annotations
#' new_con <- export_project_data(
#'   source_con, "annotated.db",
#'   tables = c("input_files", "reference_genomes", "reference_sequences", "vcf_data", 
#'              "flanking_sequences", "blast_parameters", "blast_results", "annotations")
#' )
#' }
#'
#' @export
export_project_data <- function(source_con, target_db_path, tables = NULL, force = FALSE, verbose = TRUE) {
  
  # Default table export order (respects foreign key dependencies)
  default_tables <- c(
    "input_files",
    "reference_genomes", 
    "reference_sequences",
    "vcf_data",
    "candidate_loci",
    "locus_statistics", 
    "flanking_sequences",
    "blast_parameters",
    "blast_results", 
    "annotations",
    "go_terms",
    "kegg_references",
    "pfam_domains",
    "interpro_families", 
    "eggnog_categories",
    "uniprot_cache"
  )
  
  if (is.null(tables)) tables <- default_tables
  
  if (verbose) {
    message("=== Exporting Project Data ===")
    message("Source: ", DBI::dbGetInfo(source_con)$dbname)
    message("Target: ", target_db_path)
  }
  
  # Validate source connection
  if (!DBI::dbIsValid(source_con)) {
    stop("Invalid source database connection")
  }
  
  # Create target database
  if (verbose) message("\n--- Creating target database ---")
  target_con <- create_funseq_db(target_db_path, force = force, verbose = verbose)
  
  tryCatch({
    
    # Export each table
    for (table in tables) {
      # Check if table exists in source
      if (!table %in% DBI::dbListTables(source_con)) {
        if (verbose) message("Skipping ", table, " (not found in source)")
        next
      }
      
      # Check if table has data
      count <- DBI::dbGetQuery(source_con, paste0("SELECT COUNT(*) as count FROM ", table))$count
      
      if (count > 0) {
        if (verbose) message("Exporting ", format(count, big.mark = ","), " records from ", table)
        
        # Get all data
        data <- DBI::dbGetQuery(source_con, paste0("SELECT * FROM ", table))
        
        # Bulk insert (preserving all IDs and structure)
        DBI::dbWriteTable(target_con, table, data, append = TRUE, row.names = FALSE)
      } else {
        if (verbose) message("Skipping ", table, " (empty table)")
      }
    }
    
    if (verbose) {
      message("\n=== Export completed successfully ===")
      message("New project database: ", target_db_path)
    }
    
    return(target_con)
    
  }, error = function(e) {
    # Clean up on error
    if (DBI::dbIsValid(target_con)) {
      DBI::dbDisconnect(target_con)
    }
    if (file.exists(target_db_path)) {
      file.remove(target_db_path)
      if (verbose) message("Removed incomplete database file due to error")
    }
    stop("Export failed: ", e$message)
  })
}