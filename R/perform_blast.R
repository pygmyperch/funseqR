#' Perform BLAST Search
#'
#' This function performs a BLAST search using BLAST+. It supports
#' both blastn and blastx searches against one of several specified databases.
#'
#' @param sequences A Biostrings::DNAStringSet object containing the query sequences.
#' @param db_path Character string specifying the path to the BLAST database.
#' @param db_name Character string specifying the name of the BLAST database.
#' @param blast_type Character string specifying the type of BLAST search.
#'   Must be either "blastn" or "blastx".
#' @param e_value Numeric value specifying the E-value threshold for hits.
#'   Default is 1e-5.
#' @param max_hits Integer specifying the maximum number of hits to return
#'   per query. Default is 5.
#' @param threads Integer specifying the number of CPU threads to use.
#'   Default is 1.
#' @param output_base Character string specifying the base name for output files.
#'   Default is "blast_output".
#'
#' @return A list containing two elements:
#'   \item{results}{A data frame of BLAST results, or NULL if no hits found.}
#'   \item{metadata}{A list containing metadata about the BLAST search,
#'   including parameters used, database information, and file paths.}
#'
#' @details
#' This function performs the following steps:
#' 1. Checks if BLAST+ is installed and if the specified database exists.
#' 2. Writes query sequences to a FASTA file.
#' 3. Constructs and executes the BLAST command.
#' 4. Parses the BLAST output and converts it to a data frame.
#' 5. Saves the results and metadata as a .RData file.
#'
#' The function creates several output files:
#' - A FASTA file containing the query sequences
#' - A text file containing the raw BLAST output
#' - A .RData file containing the parsed results and metadata
#'
#' @note
#' - This function requires BLAST+ to be installed and in the system PATH.
#' - The specified BLAST database must exist and have the necessary files
#'   (.nhr, .nin, .nsq for blastn; .phr, .pin, .psq for blastx).
#' - Large sequences or databases may require significant time and computational resources.
#'
#' @examples
#' \dontrun{
#' library(Biostrings)
#' sequences <- DNAStringSet(c("ATGCATGC", "GCATGCAT"))
#' blast_results <- perform_blast(
#'   sequences = sequences,
#'   db_path = "/path/to/blast/db",
#'   db_name = "nr",
#'   blast_type = "blastx",
#'   e_value = 1e-10,
#'   max_hits = 10,
#'   threads = 4,
#'   output_base = "my_blast_search"
#' )
#' }
#'
#' @importFrom Biostrings writeXStringSet
#' @importFrom utils packageVersion
#'
#' @export
perform_blast <- function(sequences, db_path, db_name,
                          blast_type = c("blastn", "blastx"),
                          e_value = 1e-5, max_hits = 5, threads = 1,
                          output_base = "blast_output") {
  blast_type <- match.arg(blast_type)

  # construct full database path
  db_file <- file.path(db_path, db_name)

  # check for .nal file (indicator of a multi-volume database)
  is_multi_volume <- file.exists(paste0(db_file, ".nal"))

  if (is_multi_volume) {
    cat("multi-volume database detected.\n")
  } else {
    # for single-volume databases, check as before
    required_extensions <- if(blast_type == "blastn") c(".nhr", ".nin", ".nsq") else c(".phr", ".pin", ".psq")
    missing_files <- sapply(required_extensions, function(ext) !file.exists(paste0(db_file, ext)))
    if (any(missing_files)) {
      stop(paste("missing required database files:",
                 paste(required_extensions[missing_files], collapse = ", "),
                 "\nplease check the path and database name."))
    }
  }

  # write sequences to fasta file
  output_fasta <- paste0(output_base, "_query_seqs.fasta")
  Biostrings::writeXStringSet(sequences, output_fasta)
  cat("query sequences written to:", output_fasta, "\n")

  # set output format
  outfmt <- "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

  # define output_blast
  output_blast <- paste0(output_base, "_", blast_type, ".txt")

  # construct blast command
  blast_command <- paste(
    blast_type,
    "-db", shQuote(db_file),
    "-query", shQuote(output_fasta),
    "-out", shQuote(output_blast),
    "-outfmt", shQuote(outfmt),
    "-evalue", e_value,
    "-max_target_seqs", max_hits,
    "-num_threads", threads
  )

  # run blast
  cat("executing blast command...\n")
  cat("command:", blast_command, "\n")
  system_result <- system(blast_command)

  if (system_result != 0) {
    stop(paste("blast command failed with exit status:", system_result))
  }

  cat("blast results written to:", output_blast, "\n")

  # read and parse blast output
  blast_output <- readLines(output_blast)
  if (length(blast_output) > 0) {
    blast_columns <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
    results <- read.table(text = blast_output, col.names = blast_columns, stringsAsFactors = FALSE)

    # convert numeric columns
    numeric_cols <- c("pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
    results[numeric_cols] <- lapply(results[numeric_cols], as.numeric)

  } else {
    warning("no blast hits found.")
    results <- NULL
  }

  # prepare the output list
  output_list <- list(
    results = results,
    command = blast_command,
    metadata = list(
      date = Sys.time(),
      blast_type = blast_type,
      db_name = db_name,
      e_value = e_value,
      max_hits = max_hits,
      threads = threads
    )
  )

  # save the output as an RDS file
  output_rds <- paste0(output_base, "_", blast_type, "_results.rds")
  saveRDS(output_list, file = output_rds)
  cat("blast results saved as RDS file:", output_rds, "\n")

  return(output_list)
}
