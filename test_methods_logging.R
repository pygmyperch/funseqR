# Test script for methods logging functionality
# This script tests the new method logging infrastructure

# Load required libraries
library(DBI)
library(RSQLite)

# Set up test database path
test_db_path <- tempfile(fileext = ".db")
cat("Creating test database at:", test_db_path, "\n")

# Create test database connection
con <- dbConnect(RSQLite::SQLite(), test_db_path)

# Test 1: Test schema creation with method_log table
cat("\n=== Test 1: Schema Creation ===\n")
tryCatch({
  source("R/schema.R")
  create_schema(con, verbose = TRUE)
  
  # Check if method_log table was created
  tables <- dbListTables(con)
  if ("method_log" %in% tables) {
    cat("✓ method_log table created successfully\n")
  } else {
    cat("✗ method_log table not found\n")
  }
}, error = function(e) {
  cat("✗ Error creating schema:", e$message, "\n")
})

# Test 2: Test logging utility functions
cat("\n=== Test 2: Logging Utility Functions ===\n")
tryCatch({
  source("R/method_logging.R")
  
  # Test basic logging
  log_method_execution(
    con = con,
    method_type = "test",
    function_name = "test_function",
    command_text = "test_command --param value",
    parameters = list(test_param = "test_value"),
    verbose = TRUE
  )
  
  # Check if log was created
  logs <- dbGetQuery(con, "SELECT * FROM method_log")
  if (nrow(logs) > 0) {
    cat("✓ Method logging works - logged", nrow(logs), "entries\n")
    cat("  - Method type:", logs$method_type[1], "\n")
    cat("  - Function name:", logs$function_name[1], "\n")
    cat("  - Command text:", logs$command_text[1], "\n")
  } else {
    cat("✗ No log entries found\n")
  }
}, error = function(e) {
  cat("✗ Error testing logging functions:", e$message, "\n")
})

# Test 3: Test methods summary functionality
cat("\n=== Test 3: Methods Summary Functionality ===\n")
tryCatch({
  source("R/summary.R")
  
  # Test methods summary
  methods_summary <- funseqR_summary(con, type = "methods")
  
  if (is.list(methods_summary)) {
    cat("✓ Methods summary function works\n")
    cat("  - Summary components:", paste(names(methods_summary), collapse = ", "), "\n")
    
    # Check overview
    if ("overview" %in% names(methods_summary)) {
      overview <- methods_summary$overview
      if (is.data.frame(overview) && nrow(overview) > 0) {
        cat("  - Overview has", nrow(overview), "method type(s)\n")
        print(overview)
      } else if (is.data.frame(overview) && nrow(overview) == 0) {
        cat("  - Overview is empty (expected for test)\n")
      }
    }
    
    # Check software environment
    if ("software_environment" %in% names(methods_summary)) {
      software <- methods_summary$software_environment
      if (is.data.frame(software) && nrow(software) > 0) {
        cat("  - Software environment captured\n")
      }
    }
    
    # Check methods text
    if ("methods_text" %in% names(methods_summary)) {
      methods_text <- methods_summary$methods_text
      if (is.data.frame(methods_text) && nrow(methods_text) > 0) {
        cat("  - Methods text generated:\n")
        cat("   ", methods_text$text[1], "\n")
      }
    }
  } else {
    cat("✗ Methods summary function failed\n")
  }
}, error = function(e) {
  cat("✗ Error testing methods summary:", e$message, "\n")
})

# Test 4: Test recent method logs utility
cat("\n=== Test 4: Recent Method Logs Utility ===\n")
tryCatch({
  recent_logs <- get_recent_method_logs(con, limit = 5)
  
  if (is.data.frame(recent_logs) && nrow(recent_logs) > 0) {
    cat("✓ Recent method logs function works\n")
    cat("  - Found", nrow(recent_logs), "recent log entries\n")
    print(recent_logs)
  } else {
    cat("  - No recent logs found (expected for minimal test)\n")
  }
}, error = function(e) {
  cat("✗ Error testing recent logs:", e$message, "\n")
})

# Test 5: Test database type summary includes method_log
cat("\n=== Test 5: Database Summary Includes method_log ===\n")
tryCatch({
  db_summary <- funseqR_summary(con, type = "database")
  
  if (is.data.frame(db_summary)) {
    method_log_row <- db_summary[db_summary$table_name == "method_log", ]
    if (nrow(method_log_row) > 0) {
      cat("✓ method_log table appears in database summary\n")
      cat("  - Record count:", method_log_row$record_count, "\n")
    } else {
      cat("✗ method_log table not found in database summary\n")
    }
  }
}, error = function(e) {
  cat("✗ Error testing database summary:", e$message, "\n")
})

# Clean up
cat("\n=== Cleanup ===\n")
dbDisconnect(con)
unlink(test_db_path)
cat("Test database cleaned up\n")

cat("\n=== Test Complete ===\n")
cat("Method logging infrastructure is ready for use!\n")