#!/bin/bash

# Clean up temporary debug and test files

echo "Cleaning up temporary debug and test files..."

# Remove debug files
rm -f debug_ora.R
rm -f debug_enrichment.R
rm -f debug_flanking_schema.R
rm -f debug_flanking_data.R
rm -f debug_export_error.R
rm -f debug_parameter_issue.R
rm -f debug_source_schema.R

# Remove check files
rm -f check_schema.R
rm -f check_flanking_schema.R

# Remove diagnose files
rm -f diagnose_flanking_issue.R

# Remove test export files
rm -f test_export_fix.R
rm -f test_export_complete_fix.R
rm -f test_debug_export.R
rm -f test_fixed_export.R
rm -f test_final_export.R
rm -f test_new_export.R
rm -f test_export_fixed.R

# Remove any test databases
rm -f test_*.db

echo "Cleanup complete!"
echo "Keeping only:"
echo "- test_schema_fix.R (final test script)"
echo "- All actual package files"