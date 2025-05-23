% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db_blast.R
\name{register_blast_params}
\alias{register_blast_params}
\title{BLAST search and results management functions for funseqR}
\usage{
register_blast_params(
  con,
  project_id,
  blast_type,
  db_name,
  db_path,
  e_value,
  max_hits,
  verbose = TRUE
)
}
\arguments{
\item{con}{A database connection object.}

\item{project_id}{The ID of the project to associate with the BLAST search.}

\item{blast_type}{Character string specifying the type of BLAST search.
Must be either "blastn" or "blastx".}

\item{db_name}{Character string specifying the name of the BLAST database.}

\item{db_path}{Character string specifying the path to the BLAST database.}

\item{e_value}{Numeric value specifying the E-value threshold for hits.}

\item{max_hits}{Integer specifying the maximum number of hits to return per query.}

\item{verbose}{Logical. If TRUE, print progress information. Default is TRUE.}
}
\value{
The ID of the newly registered BLAST parameters.
}
\description{
These functions manage BLAST parameters and results within the funseqR database.
}
\details{
Register BLAST parameters in the database

This function stores BLAST search parameters in the database.
}
