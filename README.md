# 🧬 Influenza Phylogenetic Epitope Drift Analyzer

A comprehensive bioinformatics platform for analyzing amino acid mutations at known antibody epitope sites across all 8 influenza virus segments. Designed for public health research, vaccine development, and understanding zoonotic spillover potential.

## 📋 Overview

This tool enables researchers to:

- **Track epitope drift**: Monitor amino acid changes at characterized antibody binding sites (HA1 subdomain, neuraminidase active site, etc.)
- **Identify mutation hotspots**: Find positions with high mutation rates and emerging variants
- **Build phylogenetic trees**: Visualize evolutionary relationships and sequence divergence
- **Assess conservation**: Identify stable regions suitable for vaccine targets
- **Compare your data**: Analyze custom sequences against known reference strains

## 🎯 Key Features

### 1. Multi-Segment Analysis
- Analyzes all 8 influenza segments: PB2, PB1, PA, HA, NP, NA, M, NS
- Segment-specific epitope site definitions based on literature
- Flexible selection of segments to focus on

### 2. Phylogenetic Analysis
- Constructs phylogenetic trees using UPGMA or Neighbor-Joining methods
- Calculates pairwise sequence divergence (% identity)
- Identifies temporal and geographical patterns of evolution

### 3. Epitope Tracking
- Maps mutations at known antibody binding sites:
  - **HA**: Ca, Cb, Sa, Sb epitopes; Receptor Binding Domain
  - **NA**: Active site, framework, immunodominant epitopes
  - **Other segments**: Conserved T-cell epitope regions
- Quantifies variation at each epitope position
- Highlights emerging mutations

### 4. Mutation Hotspot Detection
- Identifies positions with high mutation frequency
- Clusters nearby mutations into hotspot regions
- Distinguishes emerging variants (<30% frequency) from fixed mutations
- Calculates conservation scores for vaccine design

### 5. Temporal Drift Tracking
- Tracks mutation frequency changes over time (if collection dates available)
- Identifies acceleration patterns in epitope evolution
- Detects recently emerged variants

## 🚀 Getting Started

### Installation

```bash
# Clone repository (or extract files)
cd epitope-drift-analyzer

# Create virtual environment (recommended)
python -m venv env
source env/bin/activate  # On Windows: env\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### Running the Application

```bash
streamlit run app.py
```

The application will open at `http://localhost:8501`

## 📊 Usage Workflow

### Step 1: Load Sequences

**Option A: Upload FASTA file**
- Multi-segment FASTA format with segment names in sequence IDs
- Example format:
  ```
  >A/H3N2/HA/2023-01-15
  ATGGAGAAAATAGTGCACCCTCTCTCTC...
  >A/H3N2/NA/2023-01-15
  ATGGAGAAAATAGTGCACCCTCTCTCTC...
  ```

**Option B: NCBI GenBank IDs**
- Enter GenBank accession numbers (one per line)
- Fetches sequences automatically from NCBI using E-utilities API
- No API key required

**Option C: Sample Dataset**
- Pre-loaded H3N2 sequences spanning 2010-2024
- Useful for testing and demonstration

### Step 2: Select Analysis Segments

Choose which segments to analyze:
- **HA** (Hemagglutinin): Primary antibody target
- **NA** (Neuraminidase): Secondary antibody target
- **NP**: Internal viral protein
- **M, PB1, PB2, PA, NS**: Additional segments

### Step 3: Build Phylogenetic Trees

- Select tree construction method (UPGMA or Neighbor-Joining)
- Displays:
  - Interactive phylogenetic tree visualization
  - Pairwise sequence identity (%)
  - Maximum divergence observed
  - Alignment length and sequence count

### Step 4: Map Epitope Mutations

- Automatically identifies mutations at known antibody binding sites
- For each epitope shows:
  - Position number
  - Dominant amino acid
  - All observed variants
  - Variation percentage (0% = conserved, 100% = highly variable)
- Conservation heatmap visualization

### Step 5: Identify Mutation Hotspots

Configure hotspot detection:
- **Minimum Mutation Frequency**: Only show positions mutated in ≥N% of sequences
- **Window Size**: Clustering distance for grouping nearby mutations

Results display:
- Top mutation positions ranked by frequency
- Emerging variants (<30% frequency)
- Genomic hotspot regions
- Bar chart of mutation distribution

### Step 6: Export Results

Download analysis results:
- **CSV**: Structured data for further analysis
- **JSON**: Complete analysis including metadata
- **Report**: Text summary of findings

## 📁 Project Structure

```
epitope-drift-analyzer/
├── app.py                  # Main Streamlit application
├── phylo_analyzer.py      # Phylogenetic tree construction
├── epitope_mapper.py      # Epitope site mutation mapping
├── mutation_tracker.py    # Hotspot detection and tracking
├── ncbi_client.py         # NCBI API integration
├── requirements.txt       # Python dependencies
├── README.md             # This file
└── sample_data/          # Example sequences (if included)
```

## 🔬 Methodology

### Phylogenetic Trees
- **Multiple Sequence Alignment**: Uses MAFFT (if installed) or simple padding
- **Distance Calculation**: Percent identity between sequence pairs
- **Tree Construction**: UPGMA or Neighbor-Joining algorithm
- **Output**: Interactive Plotly visualization with branch lengths

### Epitope Mapping
- **Reference Sites**: Curated from literature for each segment
  - HA epitopes based on H3N2 crystallographic structures
  - NA epitopes from N2 active site characterization
- **Translation**: DNA sequences automatically translated to protein
- **Variation Analysis**: Position-wise amino acid frequency calculation
- **Conservation Score**: Inverse of mutation frequency

### Mutation Hotspot Detection
- **Reference Selection**: Longest sequence used as reference
- **Mutation Calling**: Position-wise comparison to reference
- **Frequency Calculation**: (Sequences with mutation / Total sequences) × 100
- **Clustering**: Groups mutations within specified window
- **Emerging Detection**: Flags variants at <30% frequency

### Temporal Drift
- **Grouping**: Sequences organized by collection date
- **Tracking**: Mutation frequency calculated per time period
- **Analysis**: Identifies acceleration or deceleration in drift

## 🧪 Example Analysis: H3N2 Drift

1. **Load**: Sample H3N2 sequences (2010-2024, 20 sequences)
2. **Select**: HA and NA segments
3. **Tree**: Shows typical antigenic drift pattern with lineage separation
4. **Epitopes**: HA Ca/Cb sites show 15-25% variation; Sa/Sb more conserved
5. **Hotspots**: Positions 137-142 (HA1) and 369-375 identified as high-variation regions
6. **Conservation**: Positions 85-95, 195-200 show >95% conservation (vaccine targets)
7. **Drift**: Clear temporal progression of mutations 2010→2024

## 📚 Citation

If using this tool for research, please cite:

```
Influenza Phylogenetic Epitope Drift Analyzer
A bioinformatics tool for tracking antibody epitope evolution
Version 1.0
```

And reference the original epitope site publications:
- **HA epitopes**: [Wiley et al., eLife 2013; Fonville et al., 2014]
- **NA epitopes**: [Yamayoshi & Kawaoka, 2019]

## 📖 Detailed Epitope Definitions

### HA (Hemagglutinin) - H3N2
- **Ca epitope** (positions 89-104): Recognized by human antibodies
- **Cb epitope** (positions 120-142): Secondary antigenic site
- **Sa epitope** (positions 195-213): Subdomain A epitope
- **Sb epitope** (positions 238-257): Subdomain B epitope
- **RBD** (positions 195-230): Receptor binding domain

### NA (Neuraminidase) - N2
- **Active site** (positions 118-129): Catalytic site region
- **Framework epitope** (positions 200-220): Structural epitope
- **Immunodominant epitope** (positions 340-360): Primary immune target

### Other Segments
- **NP**: Conserved T-cell epitopes, important for cellular immunity
- **M1**: Matrix protein, structural and immune epitopes
- **PB2/PB1/PA**: Polymerase complex, less characterized for antibody binding

## ⚙️ Advanced Usage

### Temporal Analysis
```python
from mutation_tracker import MutationTracker

tracker = MutationTracker(sequences)
temporal_df = tracker.track_temporal_drift("HA")
# Shows mutation frequency trends over time
```

### Conserved Region Identification (for vaccine design)
```python
conserved = tracker.get_conserved_regions("HA", min_conservation=95)
# Returns positions with >95% conservation across strains
```

### Emerging Variants
```python
emerging = tracker.get_emerging_mutations("HA", max_freq=20)
# Recent mutations in <20% of sequences - potential new variants
```

## 🔍 Interpreting Results

### Phylogenetic Trees
- **Branch length**: Evolutionary distance (longer = more divergence)
- **Clustering**: Sequences that group together are more similar
- **Temporal signal**: If tree shows date-based clustering, drift is directional

### Epitope Heatmap
- **Green regions**: Conserved amino acids (stable epitopes)
- **Yellow/red regions**: Variable positions (evolving epitopes)
- **Position with color intensity**: Strength of variation

### Hotspot Plots
- **Peak heights**: Mutation frequency (%) at each position
- **Clusters**: Groups of nearby mutations indicate hotspot regions
- **Color**: Conservation level (red=low conservation)

### Emerging Mutations
- **< 30% frequency**: Recently emerged, not yet fixed in population
- **Rising frequency**: Watch for potential adaptive value
- **Affects epitope**: Higher chance of escape mutant if in epitope region

## ⚠️ Limitations & Considerations

1. **Sequence Quality**: Analysis depends on input sequence accuracy
2. **Reference Dependency**: Results relative to chosen reference sequence
3. **Epitope Sites**: Based on published literature; updates needed for new discoveries
4. **Statistical Power**: Small sample sizes (<10 sequences) may show artifacts
5. **Temporal Bias**: Collection dates required for drift analysis; missing dates = no temporal analysis
6. **Alignment Assumptions**: Simple alignment; complex structural relationships may not be captured

## 🛠️ Troubleshooting

### MAFFT not found
- Fallback: Uses simple sequence padding instead
- Better alignment: Install MAFFT (`brew install mafft` or `apt-get install mafft`)

### NCBI fetch timeout
- Check internet connection
- Verify GenBank IDs are correct
- Rate limiting: waits 0.5s between requests

### Empty results
- Verify sequence format matches expected FASTA
- Check segment names in descriptions
- Try sample dataset first

## 🔐 Data Privacy

- **Local processing**: All analysis runs locally, no data sent externally (except NCBI fetch)
- **NCBI queries**: Public sequences only, no personal data submitted
- **Exports**: Saved only to your local computer

## 📝 Future Enhancements

- [ ] WebAssembly C++ alignment engine (browser-side acceleration)
- [ ] 3D protein structure visualization with epitope highlighting
- [ ] Machine learning for predicting next dominant mutations
- [ ] Integration with FluSurge for real-time surveillance
- [ ] Reassortment risk assessment (with appropriate institutional oversight)

## 📧 Support & Contribution

For questions, bug reports, or feature requests:
- Review existing documentation
- Check method parameters in source code
- Verify with sample dataset first

## 📜 License

This tool is provided for legitimate research purposes. Users are responsible for ensuring compliance with institutional biosafety and biosecurity guidelines.

---

**Built for legitimate public health phylogenetic research.**
*Last updated: 2024*
