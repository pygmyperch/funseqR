# EXPORTED

#' Database summary functions for funseqR
#'
#' These functions provide summary information about the contents of a funseqR database.
#'

#' Summarize funseqR database contents
#'
#' This function provides a summary of the data stored in a funseqR database.
#' It returns a data frame with information about each table and the number of records.
#'
#' @param con A database connection object created by \code{\link{connect_funseq_db}} or \code{\link{create_funseq_db}}.
#' @param type Character string specifying the type of summary to generate.
#'   Currently only "data" is supported, which returns record counts for all tables.
#'   Future versions will support additional types like "analyses".
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{table_name}{Character. Name of the database table.}
#'     \item{record_count}{Integer. Number of records in the table.}
#'   }
#'
#' @details
#' The function uses a template-based approach to ensure consistent output structure
#' regardless of database contents. All possible funseqR tables are included in the
#' output, with a count of 0 for tables that don't exist or are empty.
#'
#' @importFrom DBI dbListTables dbGetQuery
#'
#' @examples
#' \dontrun{
#' # Connect to a funseqR database
#' con <- connect_funseq_db("my_project.db")
#' 
#' # Get data summary
#' summary_df <- funseqR_summary(con, type = "data")
#' print(summary_df)
#' 
#' # Close connection
#' close_funseq_db(con)
#' }
#'
#' @export
funseqR_summary <- function(con, type = "data") {
  # Validate connection
  if (!DBI::dbIsValid(con)) {
    stop("Invalid database connection.")
  }
  
  # Validate type parameter
  if (!type %in% c("data")) {
    stop("Invalid type. Currently only 'data' is supported.")
  }
  
  if (type == "data") {
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
        "candidate_loci"
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
  }
  
  # Future: Add support for other summary types
  # if (type == "analyses") {
  #   # Return summary of analysis parameters, dates, etc.
  # }
}