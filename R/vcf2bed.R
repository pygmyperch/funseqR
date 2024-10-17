#' Convert loci in a VCF to BED
#'
#' This function generates a data frame in BED format and (optional) .bed file from a vcfR object.
#'
#' @title Convert VCF to BED
#' @description This function extracts the loci information from a vcfR object and returns it as a data frame. Optionally, it can also output the data as a .bed file.
#' @param vcf A vcfR object containing the VCF data.
#' @param bed_file_path Optional. The file path to the output BED file. If not supplied, the function will only return the data frame.
#' @return A data frame with columns chromosome, start_position, and end_position.
#' @importFrom vcfR getCHROM getPOS
#' @export
vcf2bed <- function(vcf, bed_file_path = NULL) {
  # Extract chromosome and position information
  chrom <- getCHROM(vcf)
  position <- getPOS(vcf)
  
  # Create the start and end positions for the .bed file
  start <- position
  end <- position + 1
  
  # Combine chromosome, start position, and end position into a data frame
  bed_data <- data.frame(chrom, start, end)
  
  # Write the data frame to a .bed file if a file path is provided
  if (!is.null(bed_file_path)) {
    write.table(bed_data, bed_file_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    message("BED file created successfully!")
  }
  
  return(bed_data)
}

# Example usage:
# vcf <- read.vcfR("path/to/your/file.vcf")
# bed_data <- vcf2bed(vcf, "path/to/output/file.bed")
# If you don't want to output to a file:
# bed_data <- vcf2bed(vcf)
