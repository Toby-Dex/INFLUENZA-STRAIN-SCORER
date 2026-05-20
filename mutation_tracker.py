"""
mutation_tracker.py
Identify mutation hotspots and emerging variants
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Tuple
from collections import defaultdict, Counter
import plotly.graph_objects as go
import plotly.express as px

class MutationTracker:
    """
    Identifies mutation hotspots and tracks emerging variants across sequences
    """
    
    def __init__(self, sequences: Dict):
        """
        Args:
            sequences: {seq_id: {sequence, segment, collection_date, ...}}
        """
        self.sequences = sequences
        self.reference = None
    
    def find_hotspots(self, segments: List[str], min_freq: float = 10, 
                     window: int = 10) -> Dict:
        """
        Identify positions with high mutation rates (mutation hotspots)
        
        Args:
            segments: segments to analyze
            min_freq: minimum mutation frequency (%) to include
            window: amino acid window for clustering mutations
        
        Returns:
            {segment: [hotspot_data]}
        """
        results = {}
        
        for segment in segments:
            segment_seqs = self._filter_by_segment(segment)
            if not segment_seqs:
                continue
            
            # Get reference sequence (longest or first)
            reference = self._get_reference(segment_seqs)
            
            # Find mutations
            mutations = self._find_mutations(segment_seqs, reference, segment)
            
            # Filter by frequency
            filtered_mutations = [m for m in mutations if m["Frequency %"] >= min_freq]
            
            # Cluster into hotspots
            hotspots = self._cluster_hotspots(filtered_mutations, window)
            
            # Sort by frequency
            hotspots.sort(key=lambda x: x["Frequency %"], reverse=True)
            
            results[segment] = hotspots
        
        return results
    
    def _filter_by_segment(self, segment: str) -> Dict:
        """Extract sequences for specific segment"""
        filtered = {}
        for seq_id, data in self.sequences.items():
            if data.get("segment") == segment:
                filtered[seq_id] = data
        return filtered
    
    def _get_reference(self, segment_seqs: Dict) -> str:
        """Select reference sequence (longest sequence)"""
        if not segment_seqs:
            return ""
        
        # Use longest sequence as reference
        longest_seq = max(segment_seqs.items(), 
                         key=lambda x: len(str(x[1].get("sequence", ""))))
        return str(longest_seq[1]["sequence"])
    
    def _find_mutations(self, segment_seqs: Dict, reference: str, 
                       segment: str) -> List[Dict]:
        """
        Find all mutations relative to reference
        Returns list of mutations with frequency
        """
        mutations = defaultdict(lambda: {"positions": [], "variants": Counter()})
        
        for seq_id, data in segment_seqs.items():
            seq = str(data.get("sequence", ""))
            
            # Align to reference (simple alignment - same length)
            if len(seq) != len(reference):
                # Pad or trim to reference length
                seq = seq[:len(reference)].ljust(len(reference), '-')
            
            # Find differences
            for pos, (ref_aa, seq_aa) in enumerate(zip(reference, seq)):
                if ref_aa != seq_aa and ref_aa != '-' and seq_aa != '-':
                    mutation_key = (pos, ref_aa, seq_aa)
                    mutations[mutation_key]["positions"].append(pos)
                    mutations[mutation_key]["variants"].update([f"{ref_aa}{pos+1}{seq_aa}"])
        
        # Convert to list and calculate frequencies
        mutation_list = []
        total_seqs = len(segment_seqs)
        
        for (pos, ref_aa, seq_aa), data in mutations.items():
            seq_with_mutation = len(data["positions"])
            freq_percent = (seq_with_mutation / total_seqs) * 100
            
            mutation_list.append({
                "Position": int(pos) + 1,  # 1-indexed
                "Ref AA": ref_aa,
                "Alt AA": seq_aa,
                "Mutation": f"{ref_aa}{pos+1}{seq_aa}",
                "Frequency %": round(freq_percent, 1),
                "Sequences": seq_with_mutation,
                "Total": total_seqs,
                "Emerging": freq_percent < 30  # Emerging if in <30% of seqs
            })
        
        # Sort by frequency
        mutation_list.sort(key=lambda x: x["Frequency %"], reverse=True)
        
        return mutation_list
    
    def _cluster_hotspots(self, mutations: List[Dict], window: int) -> List[Dict]:
        """
        Cluster nearby mutations into hotspots
        
        Args:
            mutations: list of mutation dicts
            window: amino acid distance to cluster
        
        Returns:
            clustered hotspots with aggregated frequency
        """
        if not mutations:
            return []
        
        # Sort by position
        sorted_muts = sorted(mutations, key=lambda x: x["Position"])
        
        hotspots = []
        current_cluster = []
        
        for mutation in sorted_muts:
            if not current_cluster:
                current_cluster = [mutation]
            else:
                # Check if within window of last mutation in cluster
                last_pos = current_cluster[-1]["Position"]
                if abs(mutation["Position"] - last_pos) <= window:
                    current_cluster.append(mutation)
                else:
                    # Finalize current cluster
                    hotspot = self._aggregate_cluster(current_cluster)
                    hotspots.append(hotspot)
                    current_cluster = [mutation]
        
        # Don't forget last cluster
        if current_cluster:
            hotspot = self._aggregate_cluster(current_cluster)
            hotspots.append(hotspot)
        
        return hotspots
    
    def _aggregate_cluster(self, cluster: List[Dict]) -> Dict:
        """Aggregate mutations in a cluster into hotspot data"""
        positions = [m["Position"] for m in cluster]
        
        # Average frequency
        avg_freq = np.mean([m["Frequency %"] for m in cluster])
        
        # Most common mutation
        top_mut = max(cluster, key=lambda x: x["Frequency %"])
        
        return {
            "Position": f"{min(positions)}-{max(positions)}",
            "Center Position": int(np.mean(positions)),
            "Mutations in Cluster": len(cluster),
            "Frequency %": round(avg_freq, 1),
            "Most Common": top_mut["Mutation"],
            "Emerging": any(m["Emerging"] for m in cluster),
            "Variants": ', '.join([m["Mutation"] for m in cluster[:5]])  # Top 5
        }
    
    def track_temporal_drift(self, segment: str) -> pd.DataFrame:
        """
        Track mutation frequency changes over time (if collection_date available)
        """
        segment_seqs = self._filter_by_segment(segment)
        
        # Group by collection date
        dated_seqs = defaultdict(list)
        for seq_id, data in segment_seqs.items():
            date = data.get("collection_date", "Unknown")
            dated_seqs[date].append(data)
        
        # For each date, find mutations
        temporal_data = []
        
        for date in sorted(dated_seqs.keys()):
            if date == "Unknown":
                continue
            
            seqs = dated_seqs[date]
            reference = self._get_reference({f"ref_{i}": s for i, s in enumerate(seqs)})
            
            mutations = self._find_mutations(
                {f"seq_{i}": s for i, s in enumerate(seqs)},
                reference,
                segment
            )
            
            for mut in mutations:
                temporal_data.append({
                    "Date": date,
                    "Mutation": mut["Mutation"],
                    "Frequency %": mut["Frequency %"],
                    "Sequences": mut["Sequences"]
                })
        
        return pd.DataFrame(temporal_data)
    
    def get_emerging_mutations(self, segment: str, max_freq: float = 30) -> List[Dict]:
        """
        Get recently emerged mutations (low frequency)
        
        Args:
            segment: segment to analyze
            max_freq: maximum frequency (%) to be considered "emerging"
        
        Returns:
            list of emerging mutations
        """
        segment_seqs = self._filter_by_segment(segment)
        reference = self._get_reference(segment_seqs)
        
        mutations = self._find_mutations(segment_seqs, reference, segment)
        
        # Filter for emerging
        emerging = [m for m in mutations if m["Frequency %"] <= max_freq]
        
        return emerging
    
    def get_conserved_regions(self, segment: str, min_conservation: float = 95) -> List[Dict]:
        """
        Identify highly conserved regions (few mutations)
        Useful for vaccine target identification
        
        Args:
            segment: segment to analyze
            min_conservation: minimum % conservation threshold
        
        Returns:
            list of conserved regions
        """
        segment_seqs = self._filter_by_segment(segment)
        reference = self._get_reference(segment_seqs)
        
        # Check each position
        conserved_regions = []
        window_size = 10
        total_seqs = len(segment_seqs)
        
        for start_pos in range(0, len(reference) - window_size, window_size):
            end_pos = start_pos + window_size
            window = reference[start_pos:end_pos]
            
            # Count matches in this window
            match_count = 0
            for seq_id, data in segment_seqs.items():
                seq = str(data.get("sequence", ""))
                if len(seq) >= end_pos:
                    seq_window = seq[start_pos:end_pos]
                    if seq_window == window:
                        match_count += 1
            
            conservation_pct = (match_count / total_seqs) * 100
            
            if conservation_pct >= min_conservation:
                conserved_regions.append({
                    "Position": f"{start_pos+1}-{end_pos}",
                    "Sequence": window,
                    "Conservation %": round(conservation_pct, 1),
                    "Sequences": match_count
                })
        
        return conserved_regions
