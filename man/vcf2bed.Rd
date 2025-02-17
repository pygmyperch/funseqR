% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcf2bed.R
\name{vcf2bed}
\alias{vcf2bed}
\title{Convert VCF to BED}
\usage{
vcf2bed(vcf, bed_file_path = NULL)
}
\arguments{
\item{vcf}{A vcfR object containing the VCF data.}

\item{bed_file_path}{Optional. The file path to the output BED file. If not supplied, the function will only return the data frame.}
}
\value{
A data frame with columns chromosome, start_position, and end_position.
}
\description{
This function extracts the loci information from a vcfR object and returns it as a data frame. Optionally, it can also output the data as a .bed file.
}
\details{
Convert loci in a VCF to BED

This function generates a data frame in BED format and (optional) .bed file from a vcfR object.
}
