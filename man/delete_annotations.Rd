% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db_annotation.R
\name{delete_annotations}
\alias{delete_annotations}
\title{Delete annotations from the database}
\usage{
delete_annotations(con, blast_param_id, confirm = TRUE, verbose = TRUE)
}
\arguments{
\item{con}{A database connection object.}

\item{blast_param_id}{The ID of the BLAST parameters.}

\item{confirm}{Logical. If TRUE, user confirmation will be required. Default is TRUE.}

\item{verbose}{Logical. If TRUE, print progress information. Default is TRUE.}
}
\value{
Logical. TRUE if the deletion was successful.
}
\description{
Delete annotations from the database
}
