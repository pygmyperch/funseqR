% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db_annotation.R
\name{store_uniprot_data}
\alias{store_uniprot_data}
\title{Store UniProt data in cache with robust error handling}
\usage{
store_uniprot_data(con, accession, response, debug = FALSE)
}
\arguments{
\item{con}{A database connection object}

\item{accession}{The UniProt accession number}

\item{response}{The UniProt API response}

\item{debug}{If TRUE, print debugging information}
}
\value{
TRUE if storage successful, FALSE otherwise
}
\description{
Store UniProt data in cache with robust error handling
}
