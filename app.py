"""
Sequence Analysis Platform
===========================================================================
A professional-grade dashboard for phylogenetic and epitope analysis.

Design  : Data-Dense Precision · Navy × Slate × Teal
Fonts   : JetBrains Mono (data/headings) · IBM Plex Sans (UI/body)
Charts  : Plotly with 3D visualizations where applicable
Audience: Researchers and Data Scientists
===========================================================================
"""

import json
from datetime import datetime
from collections import defaultdict, Counter
from typing import Dict, List

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import streamlit as st

# ─────────────────────────────────────────────────────────────────────────────
# PAGE CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────

st.set_page_config(
    page_title="Sequence Analysis Platform",
    page_icon=None,
    layout="wide",
    initial_sidebar_state="expanded",
)

# ─────────────────────────────────────────────────────────────────────────────
# DESIGN TOKENS
# ─────────────────────────────────────────────────────────────────────────────

PALETTE = {
    "navy":    "#0F172A",
    "navym":   "#1E3A5F",
    "blue":    "#1D4ED8",
    "teal":    "#0D9488",
    "cyan":    "#0891B2",
    "surface": "#FFFFFF",
    "bg":      "#F8FAFC",
    "text":    "#020617",
    "muted":   "#475569",
    "border":  "#CBD5E1",
    "grid":    "#E2E8F0",
    "success": "#059669",
    "warning": "#D97706",
    "danger":  "#DC2626",
}

COLORWAY = [
    "#1D4ED8",  # blue-700
    "#0D9488",  # teal-600
    "#7C3AED",  # violet-600
    "#D97706",  # amber-600
    "#DC2626",  # red-600
    "#059669",  # emerald-600
    "#0891B2",  # cyan-600
    "#64748B",  # slate-500
]

_PLOTLY_LAYOUT = dict(
    font=dict(family="IBM Plex Sans, system-ui, sans-serif", color=PALETTE["text"], size=12),
    paper_bgcolor="rgba(0,0,0,0)",
    plot_bgcolor="rgba(248,250,252,0.5)",
    colorway=COLORWAY,
    xaxis=dict(
        gridcolor=PALETTE["grid"],
        linecolor=PALETTE["border"],
        zerolinecolor=PALETTE["grid"],
        tickfont=dict(family="JetBrains Mono, monospace", size=10),
    ),
    yaxis=dict(
        gridcolor=PALETTE["grid"],
        linecolor=PALETTE["border"],
        zerolinecolor=PALETTE["grid"],
        tickfont=dict(family="JetBrains Mono, monospace", size=10),
    ),
    legend=dict(
        bgcolor="rgba(255,255,255,0.92)",
        bordercolor=PALETTE["border"],
        borderwidth=1,
        font=dict(size=11),
    ),
    hoverlabel=dict(
        bgcolor="white",
        bordercolor=PALETTE["border"],
        font=dict(family="IBM Plex Sans, sans-serif", size=12),
    ),
    margin=dict(t=48, b=44, l=52, r=24),
    title=dict(
        font=dict(family="JetBrains Mono, monospace", size=12, color=PALETTE["navy"]),
        x=0,
        xanchor="left",
        pad=dict(l=0, t=0),
    ),
)

# ─────────────────────────────────────────────────────────────────────────────
# THEME FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

def apply_theme(fig: go.Figure, title: str = "") -> go.Figure:
    """Apply professional theme to Plotly figure."""
    layout = dict(_PLOTLY_LAYOUT)
    if title:
        layout["title"] = {**layout.get("title", {}), "text": title}
    fig.update_layout(**layout)
    return fig

def _th_style() -> list:
    """Return professional table header styling."""
    return [{
        "selector": "thead th",
        "props": [
            ("background-color", PALETTE["navy"]),
            ("color", "#F1F5F9"),
            ("font-family", "JetBrains Mono, monospace"),
            ("font-size", "0.72rem"),
            ("letter-spacing", "0.06em"),
            ("padding", "8px 12px"),
        ],
    }]

# ─────────────────────────────────────────────────────────────────────────────
# STYLESHEET
# ─────────────────────────────────────────────────────────────────────────────

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=JetBrains+Mono:ital,wght@0,400;0,500;0,600;0,700;1,400&family=IBM+Plex+Sans:wght@300;400;500;600;700&display=swap');

:root {
    --navy:   #0F172A;
    --navym:  #1E3A5F;
    --blue:   #1D4ED8;
    --teal:   #0D9488;
    --bg:     #F8FAFC;
    --surf:   #FFFFFF;
    --text:   #020617;
    --muted:  #475569;
    --border: #CBD5E1;
    --grid:   #E2E8F0;
    --ok:     #059669;
    --warn:   #D97706;
    --err:    #DC2626;
}

html, body, [data-testid="stAppViewContainer"] {
    font-family: 'IBM Plex Sans', system-ui, sans-serif;
    background: var(--bg) !important;
    color: var(--text);
}

[data-testid="stSidebar"] {
    background: linear-gradient(180deg, #0F172A 0%, #1A2F4F 60%, #0C2340 100%) !important;
    border-right: 1px solid rgba(255,255,255,0.05) !important;
    box-shadow: 2px 0 16px rgba(0,0,0,0.20);
}

[data-testid="stSidebar"] * { color: #CBD5E1 !important; }

[data-testid="stSidebar"] h1,
[data-testid="stSidebar"] h2,
[data-testid="stSidebar"] h3 {
    color: #F1F5F9 !important;
    font-family: 'JetBrains Mono', monospace !important;
    font-size: 0.9rem !important;
    letter-spacing: 0.05em;
}

[data-testid="stSidebar"] label {
    color: #94A3B8 !important;
    font-size: 0.68rem !important;
    font-weight: 700 !important;
    text-transform: uppercase;
    letter-spacing: 0.10em;
}

[data-testid="stSidebar"] .stMarkdown p { 
    color: #64748B !important; 
    font-size: 0.78rem; 
}

[data-testid="stSidebar"] hr { 
    border-color: rgba(148,163,184,0.15) !important; 
}

[data-testid="stSidebar"] .stMultiSelect [data-baseweb="tag"] {
    background-color: #1D4ED8 !important;
    color: #EFF6FF !important;
    font-size: 0.72rem !important;
    border-radius: 2px !important;
}

.stTabs [data-baseweb="tab-list"] {
    gap: 2px;
    background: transparent;
    border-bottom: 1px solid var(--border);
}

.stTabs [data-baseweb="tab"] {
    background: transparent;
    border-radius: 3px 3px 0 0;
    padding: 8px 18px;
    border: 1px solid transparent;
    border-bottom: none;
    color: var(--muted) !important;
    font-family: 'IBM Plex Sans', sans-serif;
    font-size: 0.82rem;
    font-weight: 500;
    transition: color 0.15s;
    letter-spacing: 0.02em;
}

.stTabs [aria-selected="true"] {
    background: var(--surf) !important;
    color: var(--navy) !important;
    border-color: var(--border) !important;
    border-bottom-color: var(--surf) !important;
    font-weight: 600;
}

.stTabs [data-baseweb="tab"]:hover {
    color: var(--navy) !important;
    background: rgba(255,255,255,0.5) !important;
}

[data-testid="stMetric"] {
    background: var(--surf);
    border-radius: 5px;
    padding: 16px 20px !important;
    border: 1px solid var(--border);
    box-shadow: 0 1px 3px rgba(0,0,0,0.08);
}

[data-testid="stMetric"] label {
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.65rem;
    font-weight: 700;
    text-transform: uppercase;
    letter-spacing: 0.08em;
    color: var(--muted);
}

.section-heading {
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.85rem;
    font-weight: 700;
    letter-spacing: 0.08em;
    text-transform: uppercase;
    color: var(--navy);
    margin: 1.2rem 0 0.8rem 0 !important;
}

.insight-card {
    background: linear-gradient(135deg, #F0F9FF 0%, #F0FDFA 100%);
    border: 1px solid #CFFAFE;
    border-radius: 4px;
    padding: 14px 16px;
    margin: 10px 0;
    font-size: 0.88rem;
    line-height: 1.5;
    font-family: 'IBM Plex Sans', sans-serif;
}

.export-card {
    background: var(--surf);
    border: 1px solid var(--border);
    border-radius: 4px;
    padding: 16px 20px;
    margin: 12px 0;
}

.badge {
    display: inline-block;
    background: var(--blue);
    color: white;
    padding: 3px 8px;
    border-radius: 2px;
    font-size: 0.65rem;
    font-weight: 600;
    letter-spacing: 0.04em;
    margin-right: 8px;
    margin-bottom: 4px;
    font-family: 'JetBrains Mono', monospace;
}

.stDataFrame {
    font-family: 'IBM Plex Sans', sans-serif !important;
}

button {
    font-family: 'IBM Plex Sans', sans-serif !important;
    font-weight: 500 !important;
}

</style>
""", unsafe_allow_html=True)

# ─────────────────────────────────────────────────────────────────────────────
# EPITOPE DEFINITIONS
# ─────────────────────────────────────────────────────────────────────────────

EPITOPE_SITES = {
    "HA": {
        "Ca": list(range(88, 105)),
        "Cb": list(range(119, 142)),
        "Sa": list(range(194, 213)),
        "Sb": list(range(237, 257)),
    },
    "NA": {
        "active_site": list(range(117, 130)),
        "framework": list(range(199, 220)),
    },
    "NP": {
        "conserved": list(range(0, 50)),
    },
    "M": {
        "structural": list(range(0, 50)),
    }
}

SEGMENT_DESCRIPTIONS = {
    "HA": "Hemagglutinin — primary antibody target",
    "NA": "Neuraminidase — secondary antibody target",
    "NP": "Nucleoprotein — internal viral component",
    "M": "Matrix protein — structural protein",
    "PB1": "RNA polymerase subunit 1",
    "PB2": "RNA polymerase subunit 2",
    "PA": "RNA polymerase acidic protein",
    "NS": "Nonstructural proteins",
}

# ─────────────────────────────────────────────────────────────────────────────
# SESSION STATE INITIALIZATION
# ─────────────────────────────────────────────────────────────────────────────

if "sequences" not in st.session_state:
    st.session_state.sequences = {}
if "analysis_results" not in st.session_state:
    st.session_state.analysis_results = {}

# ─────────────────────────────────────────────────────────────────────────────
# UTILITY FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

def translate_dna_to_protein(dna_seq: str) -> str:
    """Translate DNA sequence to protein using standard genetic code."""
    codon_table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
    }
    
    protein = []
    dna_seq = dna_seq.upper().replace('-', '').replace('N', '')
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3]
        if len(codon) == 3:
            protein.append(codon_table.get(codon, 'X'))
    return ''.join(protein)

def parse_fasta(content: str, segments: list) -> dict:
    """Parse FASTA format content."""
    sequences = {}
    current_id = None
    current_seq = []
    
    for line in content.split('\n'):
        if line.startswith('>'):
            if current_id:
                seq = ''.join(current_seq)
                segment = identify_segment(current_id)
                if segment in segments or not segments:
                    sequences[current_id] = {
                        "sequence": seq,
                        "segment": segment,
                        "length": len(seq)
                    }
            current_id = line[1:]
            current_seq = []
        else:
            current_seq.append(line.strip())
    
    if current_id:
        seq = ''.join(current_seq)
        segment = identify_segment(current_id)
        if segment in segments or not segments:
            sequences[current_id] = {
                "sequence": seq,
                "segment": segment,
                "length": len(seq)
            }
    
    return sequences

def identify_segment(description: str) -> str:
    """Identify sequence segment from header."""
    desc_upper = description.upper()
    for segment in ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]:
        if segment in desc_upper:
            return segment
    return "Unknown"

def compute_sequence_identity(seq1: str, seq2: str) -> float:
    """Calculate percent identity between two sequences."""
    seq1 = seq1.upper().replace('-', '')
    seq2 = seq2.upper().replace('-', '')
    
    if len(seq1) != len(seq2):
        min_len = min(len(seq1), len(seq2))
        seq1 = seq1[:min_len]
        seq2 = seq2[:min_len]
    
    if len(seq1) == 0:
        return 0
    
    matches = sum(a == b for a, b in zip(seq1, seq2))
    return (matches / len(seq1)) * 100

def build_distance_matrix(sequences: Dict) -> pd.DataFrame:
    """Build pairwise distance matrix between sequences."""
    seq_ids = list(sequences.keys())
    n = len(seq_ids)
    matrix = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i, n):
            identity = compute_sequence_identity(
                sequences[seq_ids[i]]["sequence"],
                sequences[seq_ids[j]]["sequence"]
            )
            distance = 100 - identity
            matrix[i][j] = distance
            matrix[j][i] = distance
    
    return pd.DataFrame(matrix, index=seq_ids, columns=seq_ids)

def find_mutations_at_epitopes(sequences: Dict, segment: str) -> List[Dict]:
    """Identify amino acid variations at epitope sites."""
    epitopes = EPITOPE_SITES.get(segment, {})
    if not epitopes:
        return []
    
    proteins = {}
    for seq_id, data in sequences.items():
        if data.get("segment") != segment:
            continue
        
        seq = data["sequence"]
        if any(c in seq.upper() for c in "ATCG") and len(seq) % 3 == 0:
            proteins[seq_id] = translate_dna_to_protein(seq)
        else:
            proteins[seq_id] = seq.upper()
    
    if not proteins:
        return []
    
    mutations_data = []
    total_seqs = len(proteins)
    
    for epitope_name, positions in epitopes.items():
        for pos in positions:
            aa_counts = Counter()
            
            for seq_id, protein in proteins.items():
                if pos < len(protein):
                    aa = protein[pos]
                    if aa != '-' and aa != 'X':
                        aa_counts[aa] += 1
            
            if len(aa_counts) > 1:
                dominant = aa_counts.most_common(1)[0][0]
                dominant_count = aa_counts[dominant]
                variation_pct = ((total_seqs - dominant_count) / total_seqs) * 100
                
                mutations_data.append({
                    "Epitope": epitope_name,
                    "Position": int(pos) + 1,
                    "Dominant": dominant,
                    "Variation %": round(variation_pct, 1),
                    "Sequences": total_seqs
                })
    
    return sorted(mutations_data, key=lambda x: x["Variation %"], reverse=True)

def find_hotspots(sequences: Dict, segment: str, min_freq: float = 10) -> List[Dict]:
    """Identify mutation hotspots in sequences."""
    segment_seqs = {k: v for k, v in sequences.items() if v.get("segment") == segment}
    if not segment_seqs:
        return []
    
    reference = max(
        [(k, v["sequence"]) for k, v in segment_seqs.items()],
        key=lambda x: len(x[1])
    )[1]
    
    mutations = defaultdict(int)
    total_seqs = len(segment_seqs)
    
    for seq_id, data in segment_seqs.items():
        seq = data["sequence"].upper()
        ref = reference.upper()
        
        if len(seq) < len(ref):
            seq = seq + '-' * (len(ref) - len(seq))
        else:
            seq = seq[:len(ref)]
        
        for pos, (ref_aa, seq_aa) in enumerate(zip(ref, seq)):
            if ref_aa != seq_aa and ref_aa != '-' and seq_aa != '-':
                mutations[pos] += 1
    
    hotspots = []
    for pos, count in mutations.items():
        freq = (count / total_seqs) * 100
        if freq >= min_freq:
            hotspots.append({
                "Position": int(pos) + 1,
                "Frequency %": round(freq, 1),
                "Sequences": count,
                "Total": total_seqs
            })
    
    return sorted(hotspots, key=lambda x: x["Frequency %"], reverse=True)

def load_sample_data() -> Dict:
    """Load sample sequence dataset."""
    return {
        "Sample_2010_HA": {
            "sequence": "ATGGAGAAAATAGTGCACCCTCTCTCTACGACGCGGTCTGGACATGGTGGAGTGTACCGAAATCTGTAGAAAGAGTGCAGTGATCTCTATGTAGAAACACCAAACTCTGAGAATGGAATATGCTTATCCAGGAGATTTCATCGACTACGAGGAACTGAGGGAGCTATAAACAGCAGTATCATCTCAGACCAGCAGAGTGCTCAGAGCCGGCTATACTGAAAGAGCAGTCAATCTCTCTGCTACTGGGGAGCCAGATGTATGCTATCCAGGAGATTTCATCGACTACGAGGAACTGAGGGAGCTATAAACAGCAGTATCATCTCAGACCAGCAGAGTGCTCAGAGCCGGCTATACTGAAAGAGCAGTCAATCTCTCTGCTACTGGGGAGCCAGATGTATGCTATCCAGGAGATTTCATCGACTACGAGGAACTGAGGGAGCTATAAACAGCAGTATCATCTCAGACCAGCAGAGTGCTCAGAGCCGGCTATACTGAAAGAGCAGTCAATCTCTCTGCTACTGGGGAGCCAGATGTATGCTAT",
            "segment": "HA"
        },
        "Sample_2010_NA": {
            "sequence": "ATGGAGAAAATAGTGCACCCTCTCTCTACGACGCGGTCTGGACATGGTGGAGTGTACCGAAATCTGTAGAAAGAGTGCAGTGATCTCTATGTAGAAACACCAAACTCTGAGAATGGAATATGCTTATCCAGGAGATTTCATCGACTACGAGGAACTGAGGGAGCTATAAACAGCAGTATCATCTCAGACCAGCAGAGTGCTCAGAGCCGGCTATACTGAAAGAGCAGTCAATCTCTCTGCTACTGGGGAGCCAGATGTATGCTAT",
            "segment": "NA"
        },
        "Sample_2015_HA": {
            "sequence": "ATGGAGAAAATAGTGCACCCTCTCTCTACGACGCGGTCTGGACATGGTGGAGTGTACCGAAATCTGTAGAAAGAGTGCAGTGATCTCTATGTAGAAACACCAAACTCTGAGAATGGAATATGCTTATCCAGGAGATTTCATCGACTACGAGGAACTGAGGGAGCTATAAACAGCAGTATCATCTCAGACCAGCAGAGTGCTCAGAGCCGGCTATACTGAAAGAGCAGTCAATCTCTCTGCTACTGGGGAGCCAGATGTATGCTATCCAGGAGATTTCATCGACTACGAGGAACTGAGGGAGCTATAAACAGCAGTATCATCTCAGACCAGCAGAGTGCTCAGAGCCGGCTATACTGAAAGAGCAGTCAATCTCTCTGCTACTGGGGAGCCAGATGTATGCTAT",
            "segment": "HA"
        },
        "Sample_2015_NA": {
            "sequence": "ATGGAGAAAATAGTGCACCCTCTCTCTACGACGCGGTCTGGACATGGTGGAGTGTACCGAAATCTGTAGAAAGAGTGCAGTGATCTCTATGTAGAAACACCAAACTCTGAGAATGGAATATGCTTATCCAGGAGATTTCATCGACTACGAGGAACTGAGGGAGCTATAAACAGCAGTATCATCTCAGACCAGCAGAGTGCTCAGAGCCGGCTATACTGAAAGAGCAGTCAATCTCTCTGCTACTGGGGAGCCAGATGTATGCTAT",
            "segment": "NA"
        },
        "Sample_2020_HA": {
            "sequence": "ATGGAGAAAATAGTGCACCCTCTCTCTACGACGCGGTCTGGACATGGTGGAGTGTACCGAAATCTGTAGAAAGAGTGCAGTGATCTCTATGTAGAAACACCAAACTCTGAGAATGGAATATGCTTATCCAGGAGATTTCATCGACTACGAGGAACTGAGGGAGCTATAAACAGCAGTATCATCTCAGACCAGCAGAGTGCTCAGAGCCGGCTATACTGAAAGAGCAGTCAATCTCTCTGCTACTGGGGAGCCAGATGTATGCTATCCAGGAGATTTCATCGACTACGAGGAACTGAGGGAGCTATAAACAGCAGTATCATCTCAGACCAGCAGAGTGCTCAGAGCCGGCTATACTGAAAGAGCAGTCAATCTCTCTGCTACTGGGGAGCCAGATGTATGCTAT",
            "segment": "HA"
        },
        "Sample_2020_NA": {
            "sequence": "ATGGAGAAAATAGTGCACCCTCTCTCTACGACGCGGTCTGGACATGGTGGAGTGTACCGAAATCTGTAGAAAGAGTGCAGTGATCTCTATGTAGAAACACCAAACTCTGAGAATGGAATATGCTTATCCAGGAGATTTCATCGACTACGAGGAACTGAGGGAGCTATAAACAGCAGTATCATCTCAGACCAGCAGAGTGCTCAGAGCCGGCTATACTGAAAGAGCAGTCAATCTCTCTGCTACTGGGGAGCCAGATGTATGCTAT",
            "segment": "NA"
        }
    }

# ─────────────────────────────────────────────────────────────────────────────
# PAGE HEADER
# ─────────────────────────────────────────────────────────────────────────────

st.markdown("# Sequence Analysis Platform")
st.markdown(
    "Professional sequence analysis and phylogenetic characterization for research institutions"
)

# ─────────────────────────────────────────────────────────────────────────────
# SIDEBAR CONTROLS
# ─────────────────────────────────────────────────────────────────────────────

with st.sidebar:
    st.markdown("## Data Source")
    input_mode = st.radio(
        "Select input method",
        ["Upload Dataset", "Sample Dataset"],
        label_visibility="collapsed"
    )
    
    st.markdown("---")
    st.markdown("## Analysis Parameters")
    
    selected_segments = st.multiselect(
        "Segments",
        ["HA", "NA", "NP", "M", "PB1", "PB2", "PA", "NS"],
        default=["HA", "NA"],
        help="Select sequence segments for analysis"
    )
    
    min_hotspot_freq = st.slider(
        "Hotspot Frequency Threshold (%)",
        1, 50, 10,
        help="Minimum mutation frequency to identify hotspots"
    )
    
    st.markdown("---")
    st.markdown("### Quick Reference")
    st.caption(
        "HA and NA are primary antibody targets. Other segments provide structural context."
    )

# ─────────────────────────────────────────────────────────────────────────────
# TAB LAYOUT
# ─────────────────────────────────────────────────────────────────────────────

tab_dataset, tab_epitope, tab_hotspots, tab_export = st.tabs([
    "Dataset",
    "Epitope Analysis",
    "Hotspot Mapping",
    "Export"
])

# ╔═════════════════════════════════════════════════════════════════════════════╗
# ║ TAB 1 — DATASET                                                             ║
# ╚═════════════════════════════════════════════════════════════════════════════╝

with tab_dataset:
    st.markdown('<p class="section-heading">Data Input and Management</p>', unsafe_allow_html=True)
    
    col_input, col_status = st.columns([2, 1])
    
    with col_input:
        if input_mode == "Upload Dataset":
            uploaded_file = st.file_uploader(
                "Upload FASTA file",
                type=["fasta", "fa", "faa"],
                help="Multi-sequence FASTA format"
            )
            if uploaded_file:
                content = uploaded_file.getvalue().decode("utf-8")
                st.session_state.sequences = parse_fasta(content, selected_segments if selected_segments else [])
                st.success("Dataset loaded successfully")
        
        else:  # Sample Dataset
            if st.button("Load Sample Dataset", use_container_width=True):
                st.session_state.sequences = load_sample_data()
                st.success("Sample dataset loaded (6 sequences, 2010–2020)")
    
    with col_status:
        st.metric(
            "Sequences Loaded",
            len(st.session_state.sequences)
        )
    
    if st.session_state.sequences:
        st.markdown("---")
        st.markdown('<p class="section-heading">Loaded Sequences</p>', unsafe_allow_html=True)
        
        seq_summary = []
        for seq_id, data in st.session_state.sequences.items():
            seq_summary.append({
                "Identifier": seq_id,
                "Segment": data.get("segment", "Unknown"),
                "Length (bp)": data.get("length", len(data.get("sequence", "")))
            })
        
        df_summary = pd.DataFrame(seq_summary)
        st.dataframe(
            df_summary.style
            .set_table_styles(_th_style()),
            use_container_width=True,
            hide_index=True
        )
        
        st.markdown("---")
        st.markdown('<p class="section-heading">Sequence Statistics</p>', unsafe_allow_html=True)
        
        m1, m2, m3 = st.columns(3)
        with m1:
            st.metric("Total Length", f"{sum(d.get('length', 0) for d in st.session_state.sequences.values()):,} bp")
        with m2:
            segments_found = set(d.get('segment') for d in st.session_state.sequences.values())
            st.metric("Unique Segments", len(segments_found))
        with m3:
            avg_length = np.mean([d.get('length', 0) for d in st.session_state.sequences.values()])
            st.metric("Average Length", f"{int(avg_length):,} bp")

# ╔═════════════════════════════════════════════════════════════════════════════╗
# ║ TAB 2 — EPITOPE ANALYSIS                                                    ║
# ╚═════════════════════════════════════════════════════════════════════════════╝

with tab_epitope:
    if not st.session_state.sequences:
        st.info("Load a dataset in the **Dataset** tab to begin analysis.")
    else:
        st.markdown('<p class="section-heading">Epitope Site Mutations</p>', unsafe_allow_html=True)
        st.caption(
            "Tracks amino acid variations at known antibody binding sites. "
            "Higher variation percentages indicate evolving epitopes."
        )
        
        for segment in selected_segments:
            segment_seqs = {k: v for k, v in st.session_state.sequences.items() 
                          if v.get("segment") == segment}
            if not segment_seqs:
                continue
            
            with st.expander(f"**{segment}** — {SEGMENT_DESCRIPTIONS.get(segment, 'Unknown')}", 
                           expanded=(segment == "HA")):
                mutations = find_mutations_at_epitopes(segment_seqs, segment)
                
                if mutations:
                    df_mut = pd.DataFrame(mutations)
                    
                    # Metrics
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.metric("Epitope Sites", len(df_mut["Epitope"].unique()))
                    with col2:
                        max_var = df_mut["Variation %"].max()
                        st.metric("Max Variation", f"{max_var:.1f}%")
                    with col3:
                        conserved = (df_mut["Variation %"] < 10).sum()
                        st.metric("Conserved Sites", conserved)
                    
                    st.markdown("---")
                    
                    # Table
                    st.markdown("**Mutation Details**")
                    st.dataframe(
                        df_mut.style
                        .background_gradient(subset=["Variation %"], cmap="RdYlGn_r")
                        .set_table_styles(_th_style()),
                        use_container_width=True,
                        hide_index=True
                    )
                    
                    # Chart
                    st.markdown("---")
                    st.markdown("**Variation Distribution**")
                    
                    fig = px.bar(
                        df_mut.nlargest(15, "Variation %"),
                        x="Position",
                        y="Variation %",
                        color="Variation %",
                        color_continuous_scale="RdYlGn_r",
                        hover_data=["Epitope", "Dominant"],
                    )
                    apply_theme(fig, f"Epitope Site Variation — {segment}")
                    st.plotly_chart(fig, use_container_width=True)
                    
                else:
                    st.info(f"No epitope sites defined for {segment}")

# ╔═════════════════════════════════════════════════════════════════════════════╗
# ║ TAB 3 — HOTSPOT MAPPING                                                     ║
# ╚═════════════════════════════════════════════════════════════════════════════╝

with tab_hotspots:
    if not st.session_state.sequences:
        st.info("Load a dataset in the **Dataset** tab to begin analysis.")
    else:
        st.markdown('<p class="section-heading">Mutation Hotspot Detection</p>', unsafe_allow_html=True)
        st.caption(
            f"Identifying positions with ≥{min_hotspot_freq}% mutation frequency. "
            "Hotspots indicate genomic regions under selective pressure."
        )
        
        for segment in selected_segments:
            segment_seqs = {k: v for k, v in st.session_state.sequences.items() 
                          if v.get("segment") == segment}
            if not segment_seqs:
                continue
            
            with st.expander(f"**{segment}** Hotspots", expanded=(segment == "HA")):
                hotspots = find_hotspots(segment_seqs, segment, min_hotspot_freq)
                
                if hotspots:
                    df_hot = pd.DataFrame(hotspots)
                    
                    # Metrics
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.metric("Hotspots Detected", len(df_hot))
                    with col2:
                        st.metric("Peak Frequency", f"{df_hot['Frequency %'].max():.1f}%")
                    with col3:
                        avg_freq = df_hot["Frequency %"].mean()
                        st.metric("Mean Frequency", f"{avg_freq:.1f}%")
                    
                    st.markdown("---")
                    
                    # Table
                    st.markdown("**Hotspot Positions**")
                    st.dataframe(
                        df_hot.head(20).style
                        .background_gradient(subset=["Frequency %"], cmap="Reds")
                        .set_table_styles(_th_style()),
                        use_container_width=True,
                        hide_index=True
                    )
                    
                    # 3D Visualization
                    st.markdown("---")
                    st.markdown("**Spatial Distribution (3D Surface)**")
                    
                    # Create 3D surface plot
                    x_pos = df_hot["Position"].values
                    y_freq = df_hot["Frequency %"].values
                    
                    fig_3d = go.Figure(data=[
                        go.Scatter3d(
                            x=x_pos,
                            y=np.ones_like(x_pos),
                            z=y_freq,
                            mode='markers',
                            marker=dict(
                                size=8,
                                color=y_freq,
                                colorscale='Reds',
                                showscale=True,
                                colorbar=dict(title="Frequency %")
                            ),
                            text=[f"Position {int(p)}: {f:.1f}%" for p, f in zip(x_pos, y_freq)],
                            hoverinfo="text"
                        )
                    ])
                    
                    fig_3d.update_layout(
                        title=f"3D Hotspot Distribution — {segment}",
                        scene=dict(
                            xaxis_title="Sequence Position",
                            yaxis_title="",
                            zaxis_title="Frequency %",
                            bgcolor="rgba(240,245,250,0.5)"
                        ),
                        height=500
                    )
                    apply_theme(fig_3d)
                    st.plotly_chart(fig_3d, use_container_width=True)
                    
                    # Bar chart
                    st.markdown("---")
                    st.markdown("**Frequency Ranking**")
                    
                    fig = px.bar(
                        df_hot.head(15),
                        x="Position",
                        y="Frequency %",
                        color="Frequency %",
                        color_continuous_scale="Reds",
                        hover_data=["Sequences", "Total"]
                    )
                    apply_theme(fig, f"Top Hotspots — {segment}")
                    st.plotly_chart(fig, use_container_width=True)
                
                else:
                    st.info(f"No hotspots detected for {segment} at {min_hotspot_freq}% threshold")

# ╔═════════════════════════════════════════════════════════════════════════════╗
# ║ TAB 4 — EXPORT                                                              ║
# ╚═════════════════════════════════════════════════════════════════════════════╝

with tab_export:
    if not st.session_state.sequences:
        st.info("Load a dataset in the **Dataset** tab to enable export.")
    else:
        st.markdown('<p class="section-heading">Export Analysis Results</p>', unsafe_allow_html=True)
        
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Compile results
        export_data = {
            "timestamp": datetime.now().isoformat(),
            "segments_analyzed": selected_segments,
            "sequence_count": len(st.session_state.sequences),
            "epitope_analysis": {},
            "hotspot_analysis": {}
        }
        
        for segment in selected_segments:
            segment_seqs = {k: v for k, v in st.session_state.sequences.items() 
                          if v.get("segment") == segment}
            if segment_seqs:
                mutations = find_mutations_at_epitopes(segment_seqs, segment)
                hotspots = find_hotspots(segment_seqs, segment, min_hotspot_freq)
                
                export_data["epitope_analysis"][segment] = mutations
                export_data["hotspot_analysis"][segment] = hotspots
        
        # CSV Export
        st.markdown('<div class="export-card">', unsafe_allow_html=True)
        ec1, ec2 = st.columns([3, 1])
        
        with ec1:
            st.markdown("**Structured Data Export**")
            st.caption("Complete analysis results in JSON format — suitable for integration and further processing")
            st.markdown(
                '<span class="badge">JSON</span>'
                '<span class="badge">Complete Dataset</span>'
                '<span class="badge">Timestamped</span>',
                unsafe_allow_html=True,
            )
        
        with ec2:
            json_str = json.dumps(export_data, indent=2, default=str)
            st.download_button(
                "Export JSON",
                json_str,
                file_name=f"sequence_analysis_{ts}.json",
                mime="application/json",
                use_container_width=True,
            )
        
        st.markdown("</div>", unsafe_allow_html=True)
        
        # Data Preview
        st.markdown("---")
        st.markdown('<p class="section-heading">Analysis Summary</p>', unsafe_allow_html=True)
        
        summary_items = [
            f"Sequences analyzed: **{len(st.session_state.sequences)}**",
            f"Segments processed: **{len(selected_segments)}**",
            f"Analysis timestamp: **{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}**",
        ]
        
        for item in summary_items:
            st.markdown(f'<div class="insight-card">{item}</div>', unsafe_allow_html=True)

# ─────────────────────────────────────────────────────────────────────────────
# FOOTER
# ─────────────────────────────────────────────────────────────────────────────

st.markdown("---")
st.markdown(
    "<p style='text-align: center; color: #64748B; font-size: 0.75rem;'>"
    "Sequence Analysis Platform — Professional research infrastructure<br>"
    f"Generated {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
    "</p>",
    unsafe_allow_html=True
)
