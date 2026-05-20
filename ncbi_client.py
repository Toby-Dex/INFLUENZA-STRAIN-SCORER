"""
ncbi_client.py
Fetch influenza sequences from NCBI GenBank using E-utilities API
"""

import requests
import pandas as pd
from typing import Dict, List, Optional
from Bio import SeqIO
from io import StringIO
import time
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class NCBIClient:
    """
    Client for NCBI E-utilities API
    Fetches GenBank records and parses influenza sequences
    """
    
    # Base URLs for NCBI E-utilities
    EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    EFETCH_URL = f"{EUTILS_BASE}/efetch.fcgi"
    ESEARCH_URL = f"{EUTILS_BASE}/esearch.fcgi"
    
    # Influenza segment boundaries (approximate)
    SEGMENT_BOUNDARIES = {
        "PB2": (0, 2341),
        "PB1": (0, 2341),
        "PA": (0, 2233),
        "HA": (0, 1778),
        "NP": (0, 1565),
        "NA": (0, 1410),
        "M": (0, 1027),
        "NS": (0, 890),
    }
    
    def __init__(self, email: str = "research@example.com"):
        """
        Initialize NCBI client
        
        Args:
            email: Email for NCBI (required for E-utilities access)
        """
        self.email = email
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": "PhyloEpitopeAnalyzer/1.0"})
    
    def fetch_sequences(self, genbank_ids: List[str], segments: List[str]) -> Dict:
        """
        Fetch sequences from GenBank IDs
        
        Args:
            genbank_ids: list of GenBank accession numbers
            segments: list of flu segments to extract
        
        Returns:
            {seq_id: {sequence, segment, collection_date, organism, ...}}
        """
        sequences = {}
        
        for i, gb_id in enumerate(genbank_ids):
            if not gb_id.strip():
                continue
            
            try:
                logger.info(f"Fetching {gb_id} ({i+1}/{len(genbank_ids)})")
                
                # Fetch from GenBank
                record = self._fetch_record(gb_id)
                if not record:
                    continue
                
                # Extract segment information
                seq_data = self._parse_record(record, segments)
                if seq_data:
                    sequences.update(seq_data)
                
                # Rate limiting
                time.sleep(0.5)
            
            except Exception as e:
                logger.error(f"Error fetching {gb_id}: {e}")
                continue
        
        logger.info(f"Successfully fetched {len(sequences)} sequences")
        return sequences
    
    def _fetch_record(self, gb_id: str) -> Optional[object]:
        """
        Fetch single GenBank record
        
        Args:
            gb_id: GenBank accession number
        
        Returns:
            BioPython SeqRecord or None
        """
        try:
            params = {
                "db": "nucleotide",
                "id": gb_id.strip(),
                "rettype": "fasta",
                "retmode": "text",
                "tool": "phylo_analyzer",
                "email": self.email
            }
            
            response = self.session.get(self.EFETCH_URL, params=params, timeout=10)
            response.raise_for_status()
            
            # Parse FASTA response
            fasta_content = response.text
            
            # Also fetch GenBank format for metadata
            params["rettype"] = "gb"
            response_gb = self.session.get(self.EFETCH_URL, params=params, timeout=10)
            
            # Parse using BioPython
            records = list(SeqIO.parse(StringIO(fasta_content), "fasta"))
            if records:
                return records[0]
            
            return None
        
        except Exception as e:
            logger.error(f"NCBI fetch error for {gb_id}: {e}")
            return None
    
    def _parse_record(self, record, segments: List[str]) -> Dict:
        """
        Parse GenBank record to extract segment information
        
        Args:
            record: BioPython SeqRecord
            segments: list of segments to extract
        
        Returns:
            {seq_id: {sequence, segment, metadata}}
        """
        sequences = {}
        
        try:
            seq_str = str(record.seq)
            seq_length = len(seq_str)
            
            # Determine segment based on sequence length and description
            segment = self._infer_segment(record.description, seq_length)
            
            if segment not in segments:
                return sequences
            
            # Extract metadata
            metadata = self._extract_metadata(record)
            
            # Store sequence
            seq_id = f"{record.id}_{segment}"
            sequences[seq_id] = {
                "sequence": seq_str,
                "segment": segment,
                "length": seq_length,
                "accession": record.id,
                "collection_date": metadata.get("date", "Unknown"),
                "organism": metadata.get("organism", "Influenza virus"),
                "host": metadata.get("host", "Unknown"),
                "country": metadata.get("country", "Unknown"),
                "description": record.description
            }
            
            logger.info(f"Parsed {segment} ({seq_length} bp) from {record.id}")
            
        except Exception as e:
            logger.error(f"Error parsing record {record.id}: {e}")
        
        return sequences
    
    def _infer_segment(self, description: str, seq_length: int) -> str:
        """
        Infer flu segment from sequence length and description
        
        Args:
            description: GenBank record description
            seq_length: sequence length in nucleotides
        
        Returns:
            segment name (HA, NA, PB2, etc.)
        """
        desc_upper = description.upper()
        
        # Check description for segment keywords
        for segment in ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]:
            if segment in desc_upper:
                return segment
        
        # Infer from length (approximate boundaries)
        if seq_length > 2200:
            return "PB2"  # or PB1, PA
        elif 1700 < seq_length < 1800:
            return "HA"
        elif 1400 < seq_length < 1500:
            return "NA"
        elif 1500 < seq_length < 1600:
            return "NP"
        elif 800 < seq_length < 1100:
            return "M"
        else:
            return "Unknown"
    
    def _extract_metadata(self, record) -> Dict:
        """
        Extract metadata from GenBank record features
        """
        metadata = {
            "date": "Unknown",
            "organism": "Unknown",
            "host": "Unknown",
            "country": "Unknown"
        }
        
        try:
            if hasattr(record, 'features'):
                for feature in record.features:
                    if feature.type == "source":
                        qualifiers = feature.qualifiers
                        
                        # Collection date
                        if "/collection_date=" in str(feature):
                            date_val = qualifiers.get("collection_date", ["Unknown"])[0]
                            metadata["date"] = date_val
                        
                        # Organism
                        if "organism" in qualifiers:
                            metadata["organism"] = qualifiers["organism"][0]
                        
                        # Host
                        if "host" in qualifiers:
                            metadata["host"] = qualifiers["host"][0]
                        
                        # Country
                        if "country" in qualifiers:
                            metadata["country"] = qualifiers["country"][0]
        
        except Exception as e:
            logger.debug(f"Error extracting metadata: {e}")
        
        return metadata
    
    def search_influenza(self, query: str, max_results: int = 100) -> List[str]:
        """
        Search for influenza sequences
        
        Args:
            query: search query (e.g., "influenza H3N2 hemagglutinin 2023")
            max_results: maximum number of results
        
        Returns:
            list of GenBank IDs
        """
        try:
            params = {
                "db": "nucleotide",
                "term": query,
                "retmax": max_results,
                "tool": "phylo_analyzer",
                "email": self.email
            }
            
            response = self.session.get(self.ESEARCH_URL, params=params, timeout=10)
            response.raise_for_status()
            
            # Parse XML response (simplified)
            import xml.etree.ElementTree as ET
            root = ET.fromstring(response.text)
            
            ids = []
            for id_elem in root.findall(".//Id"):
                if id_elem.text:
                    ids.append(id_elem.text)
            
            logger.info(f"Found {len(ids)} sequences matching '{query}'")
            return ids
        
        except Exception as e:
            logger.error(f"Search error: {e}")
            return []
    
    def get_reference_sequences(self, subtype: str = "H3N2") -> Dict:
        """
        Get standard reference sequences for a flu subtype
        
        Args:
            subtype: flu subtype (H3N2, H1N1, etc.)
        
        Returns:
            dictionary of reference sequences for each segment
        """
        # Pre-defined reference GenBank IDs (examples)
        references = {
            "H3N2": {
                "HA": "CY113184",  # H3 HA reference
                "NA": "CY113185",  # N2 NA reference
                "PB2": "CY113180",
                "PB1": "CY113181",
                "PA": "CY113182",
                "NP": "CY113183",
                "M": "CY113186",
                "NS": "CY113187"
            },
            "H1N1": {
                "HA": "CY018079",  # H1 HA reference
                "NA": "CY018080",  # N1 NA reference
                # ... other segments
            }
        }
        
        ref_ids = references.get(subtype, {})
        
        if not ref_ids:
            logger.warning(f"No reference sequences defined for {subtype}")
            return {}
        
        # Fetch all reference sequences
        all_seqs = self.fetch_sequences(list(ref_ids.values()), list(ref_ids.keys()))
        
        return all_seqs
