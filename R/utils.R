# Internal utility functions for funseqR package
# These functions are not exported and don't generate man pages

#' Null coalescing operator
#' 
#' Returns the left-hand side if it's not NULL or NA, otherwise returns the right-hand side
#' @noRd
`%||%` <- function(x, y) {
  if (is.null(x) || (length(x) == 1 && is.na(x))) y else x
}

# Add other utility functions here as needed
# For example:
# format_number <- function(x) format(x, big.mark = ",")
# safe_division <- function(x, y) ifelse(y == 0, 0, x/y)

#' Extract longest ORF from DNA sequence
#'
#' @param dna_seq A character string or DNAString containing the DNA sequence
#' @param min_aa_length Minimum length in amino acids for a valid ORF
#' @param return_type Either "nuc" for nucleotide sequence or "aa" for amino acid sequence
#' @param verbose Logical, print debug information
#' @return Character string containing the longest ORF, or NULL if none found
#' @importFrom Biostrings DNAString translate reverseComplement
#' @noRd
extract_longest_orf <- function(dna_seq, min_aa_length = 30, return_type = "nuc", verbose = FALSE) {
  if (is.null(dna_seq) || nchar(dna_seq) == 0) {
    return(NULL)
  }
  
  # Convert to DNAString if needed
  if (is.character(dna_seq)) {
    dna_seq <- Biostrings::DNAString(dna_seq)
  }
  
  # Get reverse complement for frames 4-6
  rev_comp <- Biostrings::reverseComplement(dna_seq)
  
  longest_orf <- NULL
  longest_orf_length <- 0
  
  # Check all 6 reading frames
  sequences <- list(
    # Forward frames (1, 2, 3)
    dna_seq,
    Biostrings::subseq(dna_seq, 2),
    Biostrings::subseq(dna_seq, 3),
    # Reverse frames (4, 5, 6)
    rev_comp,
    Biostrings::subseq(rev_comp, 2),
    Biostrings::subseq(rev_comp, 3)
  )
  
  for (frame in 1:6) {
    if (length(sequences[[frame]]) < 3) next  # Skip if too short to translate
    
    tryCatch({
      # Translate to amino acids - suppress warnings for incomplete codons and fuzzy bases
      # These warnings are expected when translating raw genomic sequences
      aa_seq <- suppressWarnings(
        Biostrings::translate(sequences[[frame]], if.fuzzy.codon = "X")
      )
      aa_string <- as.character(aa_seq)
      
      # Find ORFs (sequences between start and stop codons)
      # Look for sequences starting with M and ending with *
      orfs <- gregexpr("M[^*]*", aa_string, perl = TRUE)[[1]]
      
      if (orfs[1] != -1) {  # ORFs found
        for (i in 1:length(orfs)) {
          start_pos <- orfs[i]
          orf_length <- attr(orfs, "match.length")[i]
          
          if (orf_length >= min_aa_length) {
            if (orf_length > longest_orf_length) {
              longest_orf_length <- orf_length
              
              if (return_type == "aa") {
                longest_orf <- substr(aa_string, start_pos, start_pos + orf_length - 1)
              } else {
                # Return nucleotide sequence
                nuc_start <- (start_pos - 1) * 3 + 1
                nuc_end <- nuc_start + (orf_length * 3) - 1
                longest_orf <- as.character(Biostrings::subseq(sequences[[frame]], nuc_start, nuc_end))
              }
              
              if (verbose) {
                message("Found ORF in frame ", frame, " length: ", orf_length, " aa")
              }
            }
          }
        }
      }
    }, error = function(e) {
      if (verbose) message("Error in frame ", frame, ": ", e$message)
    })
  }
  
  return(longest_orf)
}
