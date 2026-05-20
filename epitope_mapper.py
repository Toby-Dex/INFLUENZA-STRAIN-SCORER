"""
epitope_mapper.py
Map amino acid changes at known antibody epitope sites
"""

import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from typing import Dict, List, Tuple
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
import json

class EpitopeMapper:
    """
    Maps mutations at characterized antibody binding sites for influenza proteins
    """
    
    # Known epitope sites based on literature (position in mature protein)
    EPITOPE_SITES = {
        "HA": {
            "name": "Hemagglutinin",
            "sites": {
                "Ca": list(range(89, 105)),  # HA1 subdomain - Ca epitope
                "Cb": list(range(120, 142)),  # Cb epitope
                "Sa": list(range(195, 213)),  # Sa epitope
                "Sb": list(range(238, 257)),  # Sb epitope
                "RBD": list(range(195, 230)),  # Receptor binding domain
            },
            "notes": "H3N2 numbering (mature protein starts at position 1 after signal peptide removal)"
        },
        "NA": {
            "name": "Neuraminidase",
            "sites": {
                "active_site": list(range(118, 130)),  # Catalytic residues
                "framework": list(range(200, 220)),  # Framework epitope
                "immune_epitope": list(range(340, 360)),  # Immunodominant epitope
            },
            "notes": "N2 numbering"
        },
        "NP": {
            "name": "Nucleoprotein",
            "sites": {
                "conserved_region": list(range(1, 50)),
                "variable_region": list(range(200, 250)),
            },
            "notes": "Conserved T cell epitopes"
        },
        "M1": {
            "name": "Matrix Protein 1",
            "sites": {
                "core": list(range(1, 50)),
                "flexible": list(range(100, 158)),
            },
            "notes": "Structural epitopes"
        }
    }
    
    def __init__(self):
        self.codon_table = self._get_genetic_code()
    
    def map_epitopes(self, sequences: Dict, segments: List[str]) -> Dict:
        """
        Map epitope mutations for selected segments
        
        Args:
            sequences: {seq_id: {sequence, segment, ...}}
            segments: list of segments to analyze
        
        Returns:
            {segment: {mutations: [...], conservation_plot: ...}}
        """
        results = {}
        
        for segment in segments:
            segment_seqs = self._filter_by_segment(sequences, segment)
            if not segment_seqs:
                continue
            
            # Get epitope sites for this segment
            epitope_sites = self._get_epitope_sites_for_segment(segment)
            if not epitope_sites:
                results[segment] = {"mutations": [], "conservation_plot": None}
                continue
            
            # Translate DNA to protein if needed
            protein_seqs = self._translate_if_needed(segment_seqs, segment)
            
            # Analyze mutations at epitope sites
            mutations_data = self._analyze_epitope_mutations(
                protein_seqs,
                epitope_sites,
                segment
            )
            
            # Create conservation visualization
            conservation_plot = self._plot_epitope_conservation(
                mutations_data,
                segment,
                epitope_sites
            )
            
            results[segment] = {
                "mutations": mutations_data,
                "conservation_plot": conservation_plot,
                "epitope_sites": epitope_sites
            }
        
        return results
    
    def _filter_by_segment(self, sequences: Dict, segment: str) -> Dict:
        """Extract sequences for specific segment"""
        filtered = {}
        for seq_id, data in sequences.items():
            if data.get("segment") == segment:
                filtered[seq_id] = data
        return filtered
    
    def _get_epitope_sites_for_segment(self, segment: str) -> Dict:
        """Get known epitope sites for segment"""
        # Map segment names to epitope definitions
        segment_map = {
            "HA": "HA",
            "NA": "NA",
            "NP": "NP",
            "M": "M1",
            "M1": "M1"
        }
        
        epitope_key = segment_map.get(segment)
        return self.EPITOPE_SITES.get(epitope_key, {}).get("sites", {})
    
    def _translate_if_needed(self, sequences: Dict, segment: str) -> Dict:
        """
        Translate DNA sequences to protein if needed
        Assumes input is DNA sequence (ATGC)
        """
        protein_seqs = {}
        
        for seq_id, data in sequences.items():
            seq = str(data["sequence"]).upper()
            
            # Check if already protein (contains amino acid letters)
            if any(aa in seq for aa in "EFIPQZ"):
                protein_seqs[seq_id] = seq
            else:
                # Translate DNA to protein
                try:
                    # Remove any non-codon characters
                    seq = seq.replace('-', '').replace('N', '')
                    
                    # Translate
                    protein = self._translate_dna(seq)
                    protein_seqs[seq_id] = protein
                except:
                    protein_seqs[seq_id] = seq  # Keep original if translation fails
        
        return protein_seqs
    
    def _translate_dna(self, dna_seq: str) -> str:
        """Translate DNA sequence to protein"""
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
        for i in range(0, len(dna_seq) - 2, 3):
            codon = dna_seq[i:i+3]
            aa = codon_table.get(codon, 'X')
            protein.append(aa)
        
        return ''.join(protein)
    
    def _analyze_epitope_mutations(self, protein_seqs: Dict, epitope_sites: Dict, 
                                   segment: str) -> List[Dict]:
        """
        Analyze amino acid variations at epitope sites
        
        Returns list of mutations found at epitope positions
        """
        mutations_data = []
        
        for epitope_name, positions in epitope_sites.items():
            if not positions:
                continue
            
            # Collect amino acids at each position across all sequences
            position_variations = defaultdict(lambda: defaultdict(int))
            
            for seq_id, protein in protein_seqs.items():
                for pos in positions:
                    if pos < len(protein):
                        aa = protein[pos]
                        position_variations[pos][aa] += 1
            
            # Calculate variation metrics for this epitope
            for pos, aa_counts in position_variations.items():
                total_seqs = len(protein_seqs)
                
                if len(aa_counts) > 1:  # Variation exists
                    dominant_aa = max(aa_counts, key=aa_counts.get)
                    dominant_count = aa_counts[dominant_aa]
                    variation_percent = ((total_seqs - dominant_count) / total_seqs) * 100
                    
                    # Get all variants
                    variants = ', '.join([f"{aa}:{count}" for aa, count in sorted(
                        aa_counts.items(), key=lambda x: x[1], reverse=True
                    )])
                    
                    mutations_data.append({
                        "Epitope": epitope_name,
                        "Position": int(pos),
                        "Dominant AA": dominant_aa,
                        "Variants": variants,
                        "Variation %": round(variation_percent, 1),
                        "Sequences": total_seqs
                    })
        
        # Sort by variation percentage (descending)
        mutations_data.sort(key=lambda x: x["Variation %"], reverse=True)
        
        return mutations_data
    
    def _plot_epitope_conservation(self, mutations_data: List[Dict], 
                                   segment: str, epitope_sites: Dict):
        """
        Create heatmap of epitope site conservation
        """
        if not mutations_data:
            return None
        
        try:
            df = pd.DataFrame(mutations_data)
            
            # Create conservation score (100 = fully conserved, 0 = highly variable)
            df["Conservation %"] = 100 - df["Variation %"]
            
            fig = px.bar(
                df.sort_values("Position"),
                x="Position",
                y="Conservation %",
                color="Conservation %",
                color_continuous_scale="RdYlGn",
                hover_data=["Epitope", "Dominant AA", "Variants"],
                title=f"{segment} - Epitope Site Conservation",
                labels={"Conservation %": "Conservation (%)"}
            )
            
            fig.update_layout(
                height=350,
                xaxis_title="Amino Acid Position",
                yaxis_title="Conservation %",
                hovermode='x unified'
            )
            
            return fig
        
        except Exception as e:
            print(f"Plotting error: {e}")
            return None
    
    def _get_genetic_code(self) -> Dict:
        """Standard genetic code for translation"""
        return {
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
