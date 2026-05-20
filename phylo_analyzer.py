"""
phylo_analyzer.py
Phylogenetic tree construction and sequence divergence analysis
"""

import numpy as np
import pandas as pd
from Bio import SeqIO, Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
import plotly.graph_objects as go
import plotly.express as px
from io import StringIO
from collections import defaultdict
from typing import Dict, List, Tuple
import tempfile
import subprocess

class PhyloAnalyzer:
    """
    Constructs phylogenetic trees and computes divergence metrics
    """
    
    def __init__(self, sequences: Dict, method: str = "UPGMA"):
        """
        Args:
            sequences: Dict of {seq_id: {sequence, segment, collection_date}}
            method: "UPGMA" or "Neighbor-Joining"
        """
        self.sequences = sequences
        self.method = method
        self.alignments = {}
        self.distance_matrices = {}
        self.trees = {}
        
    def analyze(self, segments: List[str]) -> Dict:
        """
        Run complete phylogenetic analysis for selected segments
        
        Returns:
            Dict: {segment: {tree_plot, metrics, distance_matrix}}
        """
        results = {}
        
        for segment in segments:
            # Filter sequences by segment
            segment_seqs = self._filter_by_segment(segment)
            if not segment_seqs:
                continue
            
            # Align sequences
            alignment = self._align_sequences(segment_seqs)
            if alignment:
                self.alignments[segment] = alignment
                
                # Build tree
                tree = self._build_tree(alignment)
                if tree:
                    self.trees[segment] = tree
                
                # Calculate metrics
                metrics = self._calculate_metrics(alignment)
                
                # Create visualization
                tree_plot = self._plot_tree(tree, segment)
                
                results[segment] = {
                    "tree_plot": tree_plot,
                    "metrics": metrics,
                    "alignment_length": len(alignment[0]),
                    "seq_count": len(alignment)
                }
        
        return results
    
    def _filter_by_segment(self, segment: str) -> Dict:
        """Extract sequences for specific segment"""
        filtered = {}
        for seq_id, data in self.sequences.items():
            if data.get("segment") == segment:
                filtered[seq_id] = data
        return filtered
    
    def _align_sequences(self, segment_seqs: Dict) -> MultipleSeqAlignment:
        """
        Align sequences using MAFFT if available, else simple alignment
        """
        try:
            # Create SeqRecord objects
            records = []
            for seq_id, data in segment_seqs.items():
                seq = str(data["sequence"])
                records.append(SeqRecord(Seq(seq), id=seq_id, description=""))
            
            if len(records) < 2:
                return None
            
            # Try to use MAFFT for alignment
            try:
                return self._align_with_mafft(records)
            except:
                # Fallback: simple alignment (padding for same length)
                return self._simple_alignment(records)
        
        except Exception as e:
            print(f"Alignment error: {e}")
            return None
    
    def _align_with_mafft(self, records) -> MultipleSeqAlignment:
        """Use MAFFT for sequence alignment"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            SeqIO.write(records, f, 'fasta')
            input_file = f.name
        
        output_file = input_file.replace('.fasta', '_aligned.fasta')
        
        try:
            subprocess.run(
                ['mafft', '--quiet', input_file],
                stdout=open(output_file, 'w'),
                timeout=30
            )
            return AlignIO.read(output_file, 'fasta')
        except:
            raise Exception("MAFFT not available")
    
    def _simple_alignment(self, records) -> MultipleSeqAlignment:
        """Fallback: pad sequences to same length"""
        max_len = max(len(r.seq) for r in records)
        for record in records:
            record.seq = record.seq + Seq('-' * (max_len - len(record.seq)))
        return MultipleSeqAlignment(records)
    
    def _build_tree(self, alignment: MultipleSeqAlignment):
        """Construct phylogenetic tree"""
        try:
            calculator = DistanceCalculator('identity')
            distance_matrix = calculator.build_distance_matrix(alignment)
            
            constructor = DistanceTreeConstructor(calculator, self.method)
            tree = constructor.build_tree(alignment)
            
            return tree
        except Exception as e:
            print(f"Tree building error: {e}")
            return None
    
    def _calculate_metrics(self, alignment: MultipleSeqAlignment) -> Dict:
        """Calculate divergence and identity metrics"""
        seqs = [str(r.seq).replace('-', '') for r in alignment]
        
        # Pairwise identity
        identities = []
        for i in range(len(seqs)):
            for j in range(i+1, len(seqs)):
                identity = self._sequence_identity(seqs[i], seqs[j])
                identities.append(identity)
        
        mean_identity = np.mean(identities) if identities else 0
        max_divergence = 100 - min(identities) if identities else 0
        
        return {
            "mean_identity": round(mean_identity, 1),
            "max_divergence": round(max_divergence, 1),
            "seq_count": len(seqs),
            "align_length": len(alignment[0])
        }
    
    def _sequence_identity(self, seq1: str, seq2: str) -> float:
        """Calculate percent identity between two sequences"""
        if len(seq1) != len(seq2):
            # Align to same length
            min_len = min(len(seq1), len(seq2))
            seq1 = seq1[:min_len]
            seq2 = seq2[:min_len]
        
        matches = sum(a == b for a, b in zip(seq1, seq2))
        return (matches / len(seq1)) * 100 if len(seq1) > 0 else 0
    
    def _plot_tree(self, tree, segment: str):
        """Create Plotly visualization of phylogenetic tree"""
        try:
            # Convert tree to coordinate system for visualization
            coords = self._tree_to_coordinates(tree)
            
            fig = go.Figure()
            
            # Add branches
            for x_coords, y_coords in self._get_tree_edges(tree):
                fig.add_trace(go.Scatter(
                    x=x_coords,
                    y=y_coords,
                    mode='lines',
                    line=dict(color='#2c3e50', width=1.5),
                    hoverinfo='skip',
                    showlegend=False
                ))
            
            # Add leaf nodes (sequences)
            leaf_names = [clade.name for clade in tree.get_terminals()]
            leaf_x = [coords.get(name, (0, 0))[0] for name in leaf_names]
            leaf_y = [coords.get(name, (0, 0))[1] for name in leaf_names]
            
            fig.add_trace(go.Scatter(
                x=leaf_x,
                y=leaf_y,
                mode='markers+text',
                marker=dict(size=6, color='#e74c3c'),
                text=leaf_names,
                textposition='middle right',
                textfont=dict(size=9),
                hoverinfo='text',
                showlegend=False
            ))
            
            fig.update_layout(
                title=f"{segment} - Phylogenetic Tree",
                xaxis_title="Evolutionary Distance",
                yaxis_title="Sequences",
                height=500,
                hovermode='closest',
                margin=dict(b=20, l=5, r=200, t=40),
                plot_bgcolor='rgba(240,240,240,0.5)'
            )
            
            return fig
        
        except Exception as e:
            print(f"Plotting error: {e}")
            return None
    
    def _tree_to_coordinates(self, tree) -> Dict:
        """Convert Bio.Phylo tree to x,y coordinates for plotting"""
        coords = {}
        
        def assign_coords(clade, x=0, y=0, spread=1):
            if clade.is_terminal():
                coords[clade.name] = (x + (clade.branch_length or 0), y)
            else:
                descendents = clade.get_terminals()
                for i, desc in enumerate(descendents):
                    assign_coords(desc, x + (clade.branch_length or 0), 
                                i - len(descendents)/2, spread/2)
        
        assign_coords(tree.clade)
        return coords
    
    def _get_tree_edges(self, tree):
        """Extract edge coordinates for tree visualization"""
        edges = []
        
        def traverse(clade, x=0):
            if not clade.is_terminal():
                for child in clade.clades:
                    x_child = x + (child.branch_length or 0.1)
                    edges.append((
                        [x, x_child],
                        [child, clade]
                    ))
                    traverse(child, x_child)
        
        traverse(tree.clade)
        return edges
