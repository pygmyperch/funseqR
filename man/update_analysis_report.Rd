% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db_reports.R
\name{update_analysis_report}
\alias{update_analysis_report}
\title{Update Analysis Report}
\usage{
update_analysis_report(
  con,
  project_id,
  section = "update",
  message = NULL,
  auto_summary = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{con}{A database connection object}

\item{project_id}{The project ID}

\item{section}{The section name (e.g., "vcf_import", "blast_search", "annotation")}

\item{message}{Custom message to include}

\item{auto_summary}{Include automatic database summary. Default is TRUE}

\item{verbose}{Logical. Print progress information. Default is TRUE}
}
\description{
Adds a new section or updates existing sections in the dynamic report
}
