#' Read a VCF file
#'
#' This function reads a VCF (Variant Call Format) file into R using the vcfR package.
#'
#' @param vcf_file A character string specifying the path to the VCF file.
#' @param verbose Logical. If TRUE, progress information is printed during file reading.
#'
#' @return A vcfR object containing the VCF data.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' vcf_data <- read_vcf("path/to/your/file.vcf")
#' }
#'
#' @importFrom vcfR read.vcfR
read_vcf <- function(vcf_file, verbose = FALSE) {
  # Check if the file exists
  if (!file.exists(vcf_file)) {
    stop("The specified VCF file does not exist: ", vcf_file)
  }
  
  # Check if the file has a .vcf or .vcf.gz extension
  if (!grepl("\\.vcf(\\.gz)?$", vcf_file, ignore.case = TRUE)) {
    warning("The file does not have a .vcf or .vcf.gz extension. It may not be a valid VCF file.")
  }
  
  # Attempt to read the VCF file
  tryCatch({
    vcf_data <- vcfR::read.vcfR(vcf_file, verbose = verbose)
    return(vcf_data)
  }, error = function(e) {
    stop("Error reading VCF file: ", e$message)
  })
}
