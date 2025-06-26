#!/bin/bash

echo "=== Staging changes ==="
git add DESCRIPTION
git add NAMESPACE
git add R/clusterprofiler_enrichment.R
git add R/enrichment_comparison.R
git add examples/clusterprofiler_workflow_example.R
git add test_clusterprofiler_integration.R
git add CLUSTERPROFILER_INTEGRATION.md
git add GIT_COMMANDS.md

echo "=== Committing changes ==="
git commit -m "Add clusterProfiler integration for GO enrichment analysis

- Add run_clusterprofiler_enrichment() function for robust GO analysis
- Integrate clusterProfiler's hypergeometric test (same as original funseqR)
- Add helper functions for comparison and format conversion
- Include complete workflow example and test scripts
- Update dependencies to include clusterProfiler and DOSE
- Maintain existing workflow compatibility
- Expected to reproduce original results (29 BP terms, 2 significant at FDR < 0.1)

Key benefits:
- Uses proven statistical implementation
- Enhanced visualization capabilities
- Reduced maintenance burden
- Future-proof community support

🤖 Generated with [Claude Code](https://claude.ai/code)

Co-Authored-By: Claude <noreply@anthropic.com>"

echo "=== Showing commit status ==="
git status

echo "=== Current branch ==="
git branch --show-current