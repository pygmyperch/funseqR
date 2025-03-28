#' Query UniProt REST API
#'
#' @param uniprot_accessions character vector of uniprot accession ids
#' @param databases character vector of databases to query (interpro, go, kegg)
#' @param timeout timeout in seconds for each request
#' @param quiet if true, suppresses progress messages
#' @return list containing annotations and query information
#' @importFrom httr GET content status_code config
#' @importFrom readr read_tsv
#' @importFrom dplyr rename mutate rowwise ungroup select distinct
#' @importFrom stringr str_replace_all str_remove_all str_trim str_extract str_split
query_uniprot_api <- function(uniprot_accessions,
                              databases = c("interpro", "go", "kegg"),
                              timeout = 300,
                              quiet = FALSE) {
  # track api queries
  queries <- list()

  # validate databases
  valid_dbs <- c("interpro", "go", "kegg")
  databases <- match.arg(databases, valid_dbs, several.ok = TRUE)

  if (!quiet) {
    message(sprintf("\nprocessing %d unique uniprot accessions", length(uniprot_accessions)))
  }

  # remove duplicates and nas
  uniprot_accessions <- unique(uniprot_accessions[!is.na(uniprot_accessions)])

  # database parameters
  db_params <- list(
    interpro = list(
      database = "InterPro",
      fields = "accession,xref_interpro_full"
    ),
    go = list(
      database = "GO",
      fields = "accession,go_id,go"
    ),
    kegg = list(
      database = "KEGG",
      fields = "accession,xref_kegg"
    )
  )

  # api config
  base_url <- "https://rest.uniprot.org/uniprotkb/search"
  config <- httr::config(
    connecttimeout = timeout,
    timeout = timeout
  )

  # get basic protein info
  query <- paste0(
    "(",
    paste0("accession:", uniprot_accessions, collapse = " OR "),
    ")"
  )

  params <- list(
    query = query,
    fields = "accession,id,gene_names,ec",
    format = "tsv",
    size = 500
  )

  # make basic request
  start_time <- Sys.time()
  res <- httr::GET(
    url = base_url,
    query = params,
    config = config
  )
  end_time <- Sys.time()

  # record query details
  queries$basic <- list(
    timestamp = start_time,
    duration = as.numeric(end_time - start_time),
    url = paste0(base_url, "?", paste(names(params), params, sep = "=", collapse = "&")),
    status = httr::status_code(res),
    success = httr::status_code(res) == 200
  )

  # initialize results
  result <- list()

  # process basic response
  if (queries$basic$success) {
    content <- httr::content(res, "text", encoding = "UTF-8")
    basic_info <- readr::read_tsv(content, quote = "", show_col_types = FALSE) %>%
      rename(
        uniprot_accession = Entry,
        uniprot_id = `Entry Name`,
        gene_names = `Gene Names`,
        ec_numbers = `EC number`
      ) %>%
      mutate(
        gene_names = str_replace_all(gene_names, "\\s+", " "),  # normalize spaces
        gene_names = str_replace_all(gene_names, "\\s", "; ")   # standardize separator
      )
    result$basic_info <- basic_info
    if (!quiet) message("retrieved basic protein information")
  } else {
    message(sprintf("failed to retrieve basic protein information - status code %d",
                    queries$basic$status))
    return(list(annotations = NULL, queries = queries))
  }

  # fetch database annotations
  fetch_db_annotations <- function(db) {
    query <- paste0(
      "(",
      paste0("accession:", uniprot_accessions, collapse = " OR "),
      ") AND database:", db_params[[db]]$database
    )

    params <- list(
      query = query,
      fields = db_params[[db]]$fields,
      format = "tsv",
      size = 500
    )

    start_time <- Sys.time()
    res <- httr::GET(url = base_url, query = params, config = config)
    end_time <- Sys.time()

    queries[[db]] <<- list(
      timestamp = start_time,
      duration = as.numeric(end_time - start_time),
      url = paste0(base_url, "?", paste(names(params), params, sep = "=", collapse = "&")),
      status = httr::status_code(res),
      success = httr::status_code(res) == 200
    )

    if (queries[[db]]$success) {
      content <- httr::content(res, "text", encoding = "UTF-8")
      df <- readr::read_tsv(content, quote = "", show_col_types = FALSE)

      # process by database type
      if (db == "interpro") {
        df <- df %>%
          rename(
            uniprot_accession = Entry,
            interpro_entries = InterPro
          ) %>%
          rowwise() %>%
          mutate(
            parsed = list({
              if (is.na(interpro_entries) || interpro_entries == "") {
                list(ids = NA, descriptions = NA)
              } else {
                entry <- str_remove_all(interpro_entries, '^"|"$|\\.$')
                entries <- unlist(str_split(entry, '";"'))
                ids <- str_trim(str_extract(entries, 'IPR\\d{6}'))
                descs <- str_trim(str_remove(entries, 'IPR\\d{6};\\s*'))
                list(
                  ids = paste(ids, collapse = "; "),
                  descriptions = paste(descs, collapse = "; ")
                )
              }
            }),
            interpro_ids = parsed$ids,
            interpro_descriptions = parsed$descriptions
          ) %>%
          ungroup() %>%
          select(-interpro_entries, -parsed)

      } else if (db == "go") {
        df <- df %>%
          rename(
            uniprot_accession = Entry,
            go_ids = `Gene Ontology IDs`,
            go_terms = `Gene Ontology (GO)`
          )

      } else if (db == "kegg") {
        df <- df %>%
          rename(
            uniprot_accession = Entry,
            kegg_refs = KEGG
          )
      }

      if (!quiet) {
        n_annotated <- nrow(df)
        message(sprintf("%s: retrieved annotations for %d/%d proteins (%.1f%%)",
                        db_params[[db]]$database,
                        n_annotated,
                        length(uniprot_accessions),
                        100 * n_annotated/length(uniprot_accessions)))
      }

      return(df)
    } else {
      if (!quiet) {
        message(sprintf("request failed for %s with status code %d",
                        db_params[[db]]$database,
                        queries[[db]]$status))
      }
      return(NULL)
    }
  }

  # get annotations for each database
  for (db in databases) {
    result[[db]] <- fetch_db_annotations(db)
  }

  return(list(
    annotations = result,
    queries = queries
  ))
}
