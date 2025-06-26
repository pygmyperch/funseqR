# Git Commands to Commit clusterProfiler Integration

Run these commands from the funseqR project root directory to commit the clusterProfiler integration:

## Option 1: Run the script
```bash
chmod +x commit_cp_enrich_branch.sh
./commit_cp_enrich_branch.sh
```

## Option 2: Manual commands

### 1. Create and switch to new branch
```bash
git checkout -b cp_enrich
```

### 2. Stage all changes
```bash
git add DESCRIPTION
git add NAMESPACE  
git add R/clusterprofiler_enrichment.R
git add R/enrichment_comparison.R
git add examples/clusterprofiler_workflow_example.R
git add test_clusterprofiler_integration.R
git add CLUSTERPROFILER_INTEGRATION.md
git add commit_cp_enrich_branch.sh
git add GIT_COMMANDS.md
```

### 3. Check what will be committed
```bash
git status
```

### 4. Create commit
```bash
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
```

### 5. Push to remote repository
```bash
git push -u origin cp_enrich
```

## Files Added/Modified

### New Files:
- `R/clusterprofiler_enrichment.R` - Main clusterProfiler integration
- `R/enrichment_comparison.R` - Comparison and utility functions
- `examples/clusterprofiler_workflow_example.R` - Complete workflow example
- `test_clusterprofiler_integration.R` - Test script
- `CLUSTERPROFILER_INTEGRATION.md` - Documentation

### Modified Files:
- `DESCRIPTION` - Added clusterProfiler and DOSE dependencies
- `NAMESPACE` - Added exports for new functions

## Verification Commands

After committing, verify the changes:

```bash
# Check current branch
git branch --show-current

# View commit history
git log --oneline -5

# See what files were added
git show --name-status HEAD

# Check branch differences
git diff main..cp_enrich --name-only
```