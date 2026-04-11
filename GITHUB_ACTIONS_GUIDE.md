# GitHub Actions Quick Start Guide

## Prerequisites

1. Push all files to your GitHub repository:
   - `generate_pharmacophores.py`
   - `requirements.txt`
   - `.github/workflows/pharmacophore_generation.yml`
   - Your CSV file in `Input/` directory

2. Ensure GitHub Actions is enabled in your repository settings

## Running the Workflow

### Step 1: Navigate to Actions

1. Go to your repository on GitHub
2. Click on the **"Actions"** tab at the top
3. You should see "Generate Interaction Pharmacophores" in the workflow list

### Step 2: Trigger Workflow

1. Click on **"Generate Interaction Pharmacophores"**
2. Click the **"Run workflow"** dropdown button (top right)
3. Configure parameters (see below)
4. Click **"Run workflow"** button

### Step 3: Download Results

1. Wait for workflow to complete (green checkmark)
2. Click on the completed workflow run
3. Scroll down to **"Artifacts"** section
4. Download:
   - `pharmacophores-complete` - All JSON files
   - `processing-summary` - Logs and error report

## Configuration Examples

### Example 1: Test Run (Small CSV)

**Use case:** Testing with your 32-entry CSV

```
csv_file: Input/pdb_ligands.csv
batch_size: (leave empty)
num_batches: (leave empty)
```

**Expected runtime:** 2-5 minutes  
**Output:** 32 JSON files (some may fail if ligand not found)

### Example 2: Medium Dataset (1,000 entries)

**Use case:** Processing 1,000 PDB entries

```
csv_file: Input/pdb_ligands.csv
batch_size: 200
num_batches: 5
```

**Expected runtime:** ~5-10 minutes (parallel)  
**Output:** ~900-1000 JSON files

### Example 3: Large Dataset (200,000 entries)

**Use case:** Your full production database

```
csv_file: Input/pdb_ligands_full.csv
batch_size: 1000
num_batches: 20
```

**Expected runtime:** ~1-2 hours (20 parallel jobs)  
**Output:** ~180,000+ JSON files (some entries will fail due to missing ligands, download issues, etc.)

### Example 4: Resuming Failed Job

**Use case:** Some batches failed, rerun entire workflow

The workflow automatically **skips already-processed entries**, so you can safely rerun:

```
csv_file: Input/pdb_ligands.csv
batch_size: 1000
num_batches: 20
```

Previously generated files are preserved; only missing entries are processed.

## Understanding Results

### Successful Entry

File: `Output/1A52_EST.json`

```json
{
  "pdb_id": "1A52",
  "ligand_name": "EST",
  "num_features": 3,
  "features": [...]
}
```

**Interpretation:** Estradiol (EST) in PDB 1A52 has 3 validated interaction features

### Failed Entry

Check `error_log.json` or `error_summary.txt` in artifacts:

```json
{
  "pdb_id": "6A93",
  "ligand_name": "PEG",
  "index": 8
}
```

**Common reasons:**
- Ligand is crystallization artifact (PEG, GOL) - has no biological interactions
- Ligand name typo in CSV
- PDB download timeout

### Zero Features

File: `Output/6A93_PEG.json`

```json
{
  "pdb_id": "6A93",
  "ligand_name": "PEG",
  "num_features": 0,
  "features": []
}
```

**Interpretation:** PEG (polyethylene glycol) has no pharmacophore features - this is **correct** for non-drug molecules

## Monitoring Progress

### During Execution

1. Click on running workflow
2. Click on job name (e.g., "Generate Pharmacophores (Batch 0)")
3. Expand steps to see real-time logs
4. Look for:
   - "Processing {PDB_ID}_{Ligand}"
   - "Validated X interaction features"
   - "Successfully generated..."

### Batch Progress

Each batch shows independent progress. Example:

```
Batch 0: Processing rows 0-999 ✓ Complete
Batch 1: Processing rows 1000-1999 ⚙️ Running
Batch 2: Processing rows 2000-2999 ⚙️ Running
...
```

## Troubleshooting

### "Workflow not found"

**Solution:** Ensure `.github/workflows/pharmacophore_generation.yml` is committed to your repository

### "CSV file not found"

**Solution:** Verify CSV path in workflow input matches your repository structure

### "Conda environment creation failed"

**Solution:** This is rare. Re-run the workflow - conda servers occasionally timeout

### "All batches failed"

**Solution:**
1. Check CSV format (must have `PDB_ID,Ligand_Name` header)
2. Run test script locally first: `python test_pharmacophore_generation.py`
3. Check GitHub Actions logs for Python errors

### "Out of disk space"

**Solution:** Happens with >50k entries. Use more batches (increase `num_batches` to 50-100)

## Cost and Limits

### GitHub Actions Limits (Free tier)

- **2,000 minutes/month** for private repos
- **Unlimited** for public repos
- **20 concurrent jobs** maximum

### Estimated Usage

- **32 entries:** ~5 minutes
- **1,000 entries (20 batches):** ~20 minutes (total across all jobs)
- **200,000 entries (200 batches):** ~200 minutes

### Recommendations

For **200k entries**:
- Use **public repository** (unlimited minutes)
- Run in **20 batches of 10k each**
- Estimated total time: 2-3 hours
- Minutes used: ~200-300 (well within free tier)

## Next Steps After Generation

1. **Download artifacts**
2. **Validate**: Check `pharmacophore_index.json` for total count
3. **Backup**: Store JSON files safely
4. **Triangle decomposition**: Process JSONs to generate triplet database
5. **Search engine**: Implement KD-tree indexing on triplets

## Advanced: Manual Parallelization

For **ultra-large datasets** (>500k entries), split CSV manually:

```bash
# Split CSV into 10 files
split -l 50000 -d --additional-suffix=.csv pdb_ligands.csv chunk_

# Run workflow 10 times with different CSV files
# Input: chunk_00.csv, chunk_01.csv, ..., chunk_09.csv
```

Each workflow run processes one chunk independently.
