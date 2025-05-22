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
