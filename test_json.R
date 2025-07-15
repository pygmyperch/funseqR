#!/usr/bin/env Rscript

# Test jsonlite behavior with empty list
library(jsonlite)

cat("Testing jsonlite::toJSON with empty list:\n")

result <- toJSON(list(), auto_unbox = TRUE)
cat("Length:", length(result), "\n")
cat("Value:", result, "\n")
cat("Class:", class(result), "\n")

# Test our fix
if (length(result) == 0) {
    result <- "{}"
} else {
    result <- as.character(result)[1]
}

cat("After fix - Length:", length(result), "\n")
cat("After fix - Value:", result, "\n")
cat("After fix - Class:", class(result), "\n")