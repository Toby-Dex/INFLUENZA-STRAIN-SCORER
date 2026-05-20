# Quick Start Guide - Epitope Drift Analyzer

## Installation (5 minutes)

### 1. Install Python 3.8+
Download from https://www.python.org/downloads/ (check "Add to PATH" on Windows)

### 2. Install Dependencies
```bash
pip install -r requirements.txt
```

### 3. (Optional) Install MAFFT
For better sequence alignment quality:
- **Mac**: `brew install mafft`
- **Linux**: `apt-get install mafft`
- **Windows**: Download from https://mafft.cbrc.jp/alignment/software/

## Running the App (1 minute)

```bash
streamlit run app.py
```

Opens automatically at http://localhost:8501

## First Analysis (10 minutes)

### Option 1: Load Sample Data (Recommended)
1. Go to **Input Data** tab
2. Click "Sample Dataset" → **Load Sample H3N2 Sequences**
3. Select segments: **HA** and **NA**
4. Go to **Phylogenetics** tab
5. Click **Build Phylogenetic Trees**
6. View results in **Epitope Tracking** and **Mutation Hotspots** tabs

### Option 2: Upload Your FASTA File
1. Prepare FASTA file with segment names in descriptions:
   ```
   >sequence_id_HA
   ATGGAGAAAATAGTGCACCCTCTCTCTC...
   >sequence_id_NA
   ATGGAGAAAATAGTGCACCCTCTCTCTC...
   ```
2. **Input Data** → **Upload FASTA**
3. Select file and segments
4. Proceed as above

### Option 3: Fetch from NCBI
1. **Input Data** → **NCBI GenBank IDs**
2. Enter GenBank accession numbers (one per line):
   ```
   CY073894
   CY073895
   CY073896
   ```
3. Click **Fetch from NCBI**
4. Proceed as above

## Understanding Results

### Phylogenetic Trees
- **Y-axis**: Individual sequences
- **X-axis**: Evolutionary distance
- **Branch length**: How different sequences are
- **Clustering**: Similar sequences group together

### Epitope Mutation Table
- **Epitope**: Which antibody binding site
- **Position**: Location in protein
- **Variation %**: % sequences with mutation (0% = conserved)
- **Red rows**: High variation (evolving epitopes)
- **Green rows**: Conserved sites (stable targets)

### Hotspot Bar Chart
- **Height**: Mutation frequency at position
- **Position**: Location in sequence
- **Emerging**: Recently evolved variants (<30% frequency)

## Common Tasks

### Find mutations at specific epitope
- Look at **Epitope Tracking** table
- Sort by "Variation %" (highest first)
- Red highlights show high variation

### Find emerging mutations (new variants)
- Go to **Mutation Hotspots**
- Set "Min. Mutation Frequency" to 5-20%
- Look for positions marked "Emerging"

### Save results
- Go to **Export Results**
- Download as CSV (Excel), JSON, or Report

### Compare custom sequences to reference
- **Input Data** → Upload FASTA
- App automatically detects and analyzes all segments

## Tips & Tricks

### Speed up analysis
- Start with 1-2 segments (HA + NA) instead of all 8
- Use "Sample Dataset" for testing

### Better alignment
- Install MAFFT: gives more accurate trees
- Ensure sequences are DNA (not protein)

### Understand epitope positions
- HA epitopes focus on HA1 subdomain (positions 1-330)
- NA epitopes focus on active site and framework

### Export for publication
- Use CSV export for tables
- Use JSON for complete computational reproducibility

## Troubleshooting

| Problem | Solution |
|---------|----------|
| "MAFFT not found" | Install MAFFT (see above) or use sample dataset |
| "NCBI fetch timeout" | Check internet, verify GenBank IDs, try later |
| "No results" | Check segment names match exactly (HA, NA, etc.) |
| "Empty trees" | Try sample dataset first; verify FASTA format |

## Next Steps

1. **Try sample data first** to understand workflow
2. **Upload your own sequences** (FASTA or GenBank IDs)
3. **Export results** for downstream analysis
4. **Read full README.md** for detailed methodology

## Getting Help

- Check **README.md** for detailed documentation
- Review source code comments in `*.py` files
- Verify your FASTA format matches example
- Test with sample dataset if having issues

---

**Ready to analyze? Run `streamlit run app.py` now!**
