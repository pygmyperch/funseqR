#' Extract flanking sequences for SNPs
#'
#' This function extracts flanking sequences for SNPs from a reference genome.
#' SNP positions can be provided either as a .bed file or a data frame.
#' If a SNP is closer to the chromosome end than the specified flank size,
#' the flanking sequence will be truncated to fit within the chromosome.
#'
#' @param snp_positions Either a character string specifying the path to a .bed file,
#'                      or a data frame with at least 3 columns: chromosome, start, and end positions.
#' @param reference_genome A DNAStringSet object containing the reference genome sequences.
#' @param flank_size The size of the flanking region to extract on each side of the SNP (default: 10000).
#'
#' @return A DNAStringSet object containing the flanking sequences for each SNP.
#'
#' @export
#'
#' @importFrom Biostrings DNAStringSet subseq
#' @importFrom utils read.table
extract_flanking_sequences <- function(snp_positions, reference_genome, flank_size = 10000) {
  if (!inherits(reference_genome, "DNAStringSet")) {
    stop("reference_genome must be a DNAStringSet object")
  }
  
  if (is.character(snp_positions) && length(snp_positions) == 1) {
    if (!file.exists(snp_positions)) {
      stop("The specified .bed file does not exist: ", snp_positions)
    }
    bed_data <- read.table(snp_positions, header = FALSE, stringsAsFactors = FALSE)
  } else if (is.data.frame(snp_positions)) {
    bed_data <- snp_positions
  } else {
    stop("snp_positions must be either a file path or a data frame")
  }
  
  if (ncol(bed_data) < 3) {
    stop("The SNP positions data should have at least 3 columns: chromosome, start, end")
  }
  
  bed_data <- bed_data[, 1:3]
  colnames(bed_data) <- c("chrom", "start", "end")
  bed_data$start <- as.numeric(bed_data$start)
  bed_data$end <- as.numeric(bed_data$end)
  
  # Extract the first part of the chromosome names from the reference genome
  ref_chrom_names <- sapply(strsplit(names(reference_genome), "\\s+"), `[`, 1)
  
  flanking_seqs <- vector("list", nrow(bed_data))
  valid_indices <- logical(nrow(bed_data))
  
  for (i in 1:nrow(bed_data)) {
    chrom <- bed_data$chrom[i]
    start <- bed_data$start[i]
    end <- bed_data$end[i]
    
    # Find the matching chromosome in the reference genome
    chrom_index <- which(ref_chrom_names == chrom)
    
    if (length(chrom_index) == 0) {
      warning(paste("Chromosome", chrom, "not found in reference genome. Skipping this SNP."))
      next
    }
    
    chrom_seq <- reference_genome[[chrom_index]]
    chrom_length <- length(chrom_seq)
    
    if (start > chrom_length || end > chrom_length || start < 1 || end < 1) {
      warning(paste("SNP position", start, "-", end, "is outside the range of chromosome", chrom, ". Skipping this SNP."))
      next
    }
    
    flank_start <- max(1, start - flank_size)
    flank_end <- min(chrom_length, end + flank_size)
    
    tryCatch({
      flanking_seq <- subseq(chrom_seq, start = flank_start, end = flank_end)
      flanking_seqs[[i]] <- flanking_seq
      valid_indices[i] <- TRUE
    }, error = function(e) {
      warning(paste("Error extracting flanking sequence for SNP at", chrom, ":", start, "-", end, ". Skipping this SNP."))
    })
  }
  
  valid_seqs <- flanking_seqs[valid_indices]
  
  if (length(valid_seqs) == 0) {
    stop("No valid flanking sequences could be extracted.")
  }
  
  result <- DNAStringSet(valid_seqs)
  
  names(result) <- paste0(
    bed_data$chrom[valid_indices], ":",
    bed_data$start[valid_indices], "-",
    bed_data$end[valid_indices]
  )
  
  return(result)
}
