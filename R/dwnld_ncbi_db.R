#' Generate scripts to download BLAST databases
#'
#' This function generates scripts to download BLAST databases
#' using various methods (wget, curl, rsync, or lftp).
#'
#' @param db_type Character string specifying the database type.
#'   Options are "nt", "nr", or "swissprot".
#' @param method Character string specifying the download method.
#'   Options are "wget", "curl", "rsync", or "lftp".
#' @param pget_n Integer specifying the number of simultaneous connections
#'   per file for lftp method. Default is 4.
#' @param parallel_downloads Integer specifying the number of files to
#'   download in parallel. Default is 4.
#' @param output_script Character string specifying the filename for the
#'   output shell script. If NULL, no script is generated but commands to 
#'   download each file are returned.
#' @param debug Logical indicating whether to print debug information.
#'   Default is FALSE.
#'
#' @return A list containing metadata about the database, file information,
#'   and download commands.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate download commands for the 'nr' database using lftp
#' nr_download <- dwnld_ncbi_db("nr", method = "lftp", 
#'                                  parallel_downloads = 4, pget_n = 8,
#'                                  output_script = "download_nr.sh")
#'
#' # View the generated commands
#' print(nr_download$download_commands$lftp)
#'
#' # Make the script executable
#' system("chmod +x download_nr.sh")
#'
#' # Run the generated script
#' system("./download_nr.sh")
#' }
dwnld_ncbi_db <- function(db_type = c("nt", "nr", "swissprot"), 
                              method = c("wget", "curl", "rsync", "lftp"),
                              pget_n = 4,
                              parallel_downloads = 4,
                              output_script = NULL,
                              debug = FALSE) {
  # match arguments
  db_type <- match.arg(db_type)
  method <- match.arg(method)
  
  # select appropriate json url
  url <- switch(db_type,
                "nt" = "https://ftp.ncbi.nlm.nih.gov/blast/db/nt-nucl-metadata.json",
                "nr" = "https://ftp.ncbi.nlm.nih.gov/blast/db/nr-prot-metadata.json",
                "swissprot" = "https://ftp.ncbi.nlm.nih.gov/blast/db/v5/swissprot-prot-metadata.json")
  
  # fetch and parse json data
  response <- GET(url)
  
  if (status_code(response) == 200) {
    json_data <- rawToChar(response$content)
    parsed_data <- fromJSON(json_data, flatten = TRUE)
    
    if(debug) {
      cat("parsed data structure:\n")
      print(str(parsed_data))
    }
    
    # extract file info and construct commands
    db_files <- parsed_data$files
    file_names <- basename(db_files)
    ftp_base_url <- "ftp://ftp.ncbi.nlm.nih.gov/blast/db/"
    rsync_base_url <- "rsync://ftp.ncbi.nlm.nih.gov/blast/db/"
    
    wget_commands <- paste("wget", db_files)
    curl_commands <- paste("curl -O", db_files)
    rsync_commands <- paste("rsync -avzP", paste0(rsync_base_url, file_names), ".")
    
    # construct a single lftp command for all files
    lftp_commands <- c(
      "lftp -c '",
      paste0("open ", ftp_base_url, ";"),
      "set cmd:parallel 1;",  # ensure we don't limit parallelism here
      paste0("set net:connection-limit ", parallel_downloads * pget_n, ";"),
      paste0("set net:max-retries 10;"),
      paste0("set net:reconnect-interval-base 5;"),
      paste0("set net:reconnect-interval-max 60;"),
      paste0("mirror -c --use-pget-n=", pget_n),
      paste0("--parallel=", parallel_downloads),
      paste("--include-glob", paste0("'", db_type, ".*.tar.gz'"), "-v"),
      "'"
    )
    lftp_commands <- paste(lftp_commands, collapse = " ")
    
    # prepare output
    result <- list(
      metadata = list(
        dbname = parsed_data$dbname,
        version = parsed_data$version,
        description = parsed_data$description,
        last_updated = parsed_data$`last-updated`,
        total_size_bytes = parsed_data$`bytes-total`,
        total_size_gb = round(parsed_data$`bytes-total` / 1e9, 2),
        number_of_sequences = parsed_data$`number-of-sequences`,
        number_of_volumes = parsed_data$`number-of-volumes`
      ),
      files = data.frame(
        name = file_names,
        url = db_files,
        stringsAsFactors = FALSE
      ),
      download_commands = list(
        wget = wget_commands,
        curl = curl_commands,
        rsync = rsync_commands,
        lftp = lftp_commands
      )
    )
    
    # add taxdb download commands
    result$taxdb_commands <- list(
      wget = "wget https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz",
      curl = "curl -O https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz",
      rsync = "rsync -avzP rsync://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz .",
      lftp = paste0("lftp -c 'open ", ftp_base_url, "; pget -n ", pget_n, " taxdb.tar.gz'")
    )
    
    # generate shell script if requested
    if (!is.null(output_script)) {
      script_content <- c(
        "#!/bin/bash",
        "",
        "# create and change to database directory",
        paste0("mkdir -p ", db_type),
        paste0("cd ", db_type),
        "",
        if (method != "lftp") {
          c(
            "# function to download files in parallel",
            "download_parallel() {",
            paste0("  parallel_downloads=", parallel_downloads),
            "  for cmd in \"$@\"; do",
            "    while [ $(jobs -p | wc -l) -ge $parallel_downloads ]; do",
            "      sleep 1",
            "    done",
            "    eval $cmd &",
            "  done",
            "  wait",
            "}",
            "",
            "# download blast database files",
            "download_parallel \\",
            paste0("  ", result$download_commands[[method]], " \\")
          )
        } else {
          c(
            "# download blast database files",
            result$download_commands$lftp
          )
        },
        "",
        "# download taxonomy database",
        "cd ..",  # move back to parent directory
        result$taxdb_commands[[method]],
        "",
        paste0("echo \"", db_type, " database download completed.\"")
      )
      writeLines(script_content, output_script)
      message("shell script generated: ", output_script)
    }
    
    return(result)
  } else {
    stop("failed to fetch the json. status code: ", status_code(response))
  }
}
