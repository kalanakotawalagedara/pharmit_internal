# Pharmit Database Creation - Quick Start

## Local Testing (2 Molecules)

### 1. Generate Pharmacophore JSONs
```bash
python generate_pharmacophores.py --csv Input/pdb_ligands.csv --output Output
```

**Expected Output:**
```
Loaded 2 entries from Input/pdb_ligands.csv
Processing 1A52_EST... ✓ (6 features)
Processing 2HYY_STI... ✓ (8 features)
Success: 2, Failed: 0
```

### 2. Build Database
```bash
python build_pharmacophore_database.py --json-dir Output --output Database --log-level INFO
```

**Expected Output:**
```
PHARMACOPHORE DATABASE BUILDER
======================================================================
Input directory:  Output
Output directory: Database
Database:         Database\pharmacophores.db

Found 2 JSON files
Processing 1A52_EST (6 features) ✓
Processing 2HYY_STI (8 features) ✓

======================================================================
DATABASE STATISTICS
======================================================================
Molecules processed:    2
Total triplets:         76
Total vectors:          104
Avg features/molecule:  7.0
Avg triplets/molecule:  38.0

Top triplet types:
  (HydrogenDonor, Hydrophobic, Hydrophobic): 24
  (HydrogenDonor, HydrogenAcceptor, Hydrophobic): 16
  (HydrogenAcceptor, Hydrophobic, Hydrophobic): 12
```

### 3. Verify Database
```bash
python verify_database.py
```

**Expected Output:**
```
DATABASE STRUCTURE
======================================================================
Tables (5):
  - molecules
  - pharma_types  
  - triplets
  - vectors
  - sqlite_sequence

Row counts:
  molecules: 2
  pharma_types: 6
  triplets: 76
  vectors: 104

✓ Database verification complete!
```

---

## Production Deployment (200k Molecules)

### 1. Prepare Large CSV
Create `Input/pdb_ligands_200k.csv` with 200,000 entries:
```csv
PDB_ID,Ligand_Name
1A52,EST
2HYY,STI
...
```

### 2. Update and Commit
```bash
git add .
git commit -m "Add 200k molecule database build"
git push origin main
```

### 3. Run GitHub Actions Workflow
1. Go to **Actions** tab in GitHub
2. Select **"Pharmacophore Generation"** workflow
3. Click **"Run workflow"**
4. Enter parameters:
   - CSV file path: `Input/pdb_ligands_200k.csv`
   - Number of batches: `10`
   - Batch size: `20000`

### 4. Monitor Progress
Watch workflow execution:
- **Setup:** Calculates batch parameters (~30 seconds)
- **Generate:** Processes 10 parallel batches (~10-20 hours)
- **Merge:** Combines all JSONs (~5 minutes)
- **Build-Database:** Creates SQLite database (~30-60 minutes)

### 5. Download Results
After workflow completes (~12-21 hours):
- Click on completed workflow run
- Download artifact: **"pharmacophore-database"**
- Extract files:
  - `pharmacophores.db` (~500 MB)
  - `molecule_index.csv` (PDB lookup table)
  - `build_log.json` (statistics)

---

## Output Files Explained

### pharmacophores.db
SQLite database with:
- **molecules:** PDB ID, ligand name, feature count
- **triplets:** Triangle decomposition with distances/angles
- **vectors:** 3D coordinates for each triplet vertex
- **pharma_types:** Type ID to name mapping

**Usage:** Pharmacophore similarity search (future implementation)

### molecule_index.csv
Fast PDB lookup table:
```csv
mol_id,pdb_id,ligand_name,chain,num_features,json_path
1,1A52,EST,A,6,Output/1A52_EST.json
2,2HYY,STI,A,8,Output/2HYY_STI.json
```

**Usage:** Map search results to PDB entries

### build_log.json
Build statistics:
```json
{
  "build_date": "2026-04-12T18:09:20.442Z",
  "database_path": "Database/pharmacophores.db",
  "statistics": {
    "molecules_processed": 200000,
    "molecules_failed": 152,
    "molecules_skipped": 0,
    "triplets_generated": 7000000,
    "vectors_stored": 21000000,
    "errors": [...]
  }
}
```

**Usage:** Quality control, error analysis

---

## Command Reference

### Generate Pharmacophores
```bash
# Full dataset
python generate_pharmacophores.py --csv <CSV_FILE> --output Output

# Batched processing (for GitHub Actions)
python generate_pharmacophores.py \
  --csv <CSV_FILE> \
  --output Output \
  --batch-start 0 \
  --batch-size 20000
```

### Build Database
```bash
# Standard build
python build_pharmacophore_database.py \
  --json-dir Output \
  --output Database

# Custom paths and logging
python build_pharmacophore_database.py \
  --json-dir Output \
  --output Database \
  --database-name pharmacophores.db \
  --index-name molecule_index.csv \
  --log-name build_log.json \
  --log-level INFO
```

### Verify Database
```bash
# Full verification
python verify_database.py

# Query specific data
python -c "import sqlite3; \
  conn = sqlite3.connect('Database/pharmacophores.db'); \
  print(conn.execute('SELECT COUNT(*) FROM molecules').fetchone()[0])"
```

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| `ImportError: No module named 'rdkit'` | Run `conda install -c conda-forge rdkit` |
| `FileNotFoundError: Input/pdb_ligands.csv` | Create CSV file or use correct path |
| Database already exists warning | Expected - duplicates are skipped automatically |
| JSON decode errors in log | Some PDBs unavailable - check `pharmacophore_generation.log` |
| GitHub Actions timeout | Reduce batch size or increase number of batches |
| Memory error during indexing | Run database build locally with more RAM |

---

## Expected Performance

### Local (2 Molecules)
- **JSON Generation:** ~5 seconds
- **Database Build:** ~1 second
- **Total:** ~6 seconds
- **Database Size:** ~50 KB

### Production (200k Molecules)
- **JSON Generation:** ~10-20 hours (parallel)
- **Database Build:** ~30-60 minutes (sequential)
- **Total:** ~12-21 hours
- **Database Size:** ~500 MB

---

## Next Steps

1. ✅ **Phase 1:** JSON generation (Complete)
2. ✅ **Phase 2:** Database building (Complete)
3. ⏳ **Phase 3:** Search implementation (Future work)
   - Query parser
   - Triangle matching
   - Recursive backtracking
   - RMSD calculation

---

**Ready to run!** Start with local testing, then scale to production.
