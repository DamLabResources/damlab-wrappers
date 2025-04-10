"""Wrapper for clipping and orienting sequences using minimap2.

This wrapper provides a Snakemake interface for the ClipQS tool, which aligns sequences
to a reference, clips unaligned regions, and ensures proper orientation.
"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2024"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.0.0"

import os
import sys
from typing import List, Optional, Tuple
import mappy as mp # type: ignore

if "snakemake" not in locals():
    import snakemake  # type: ignore

# Check version compatibility
required_version = snakemake.params.get("version", __version__)
if required_version != __version__:
    print(f"Warning: Wrapper version {__version__} does not match required version {required_version}")

def get_alignment_bounds(hits: List[mp.Alignment]) -> Optional[Tuple[int, int, float, int]]:
    """
    Get the outermost alignment bounds from all hits.
    
    Args:
        hits: List of alignment hits from minimap2
        
    Returns:
        Tuple of (start, end, coverage, strand) or None if no valid alignments
    """
    if not hits:
        return None
        
    # Get total sequence length from first hit
    seq_len = hits[0].ctg_len
    
    # Filter for primary alignments and sort by query start
    primary_hits = [h for h in hits if h.is_primary]
    if not primary_hits:
        primary_hits = hits  # Use all hits if no primary alignments
    
    hits = sorted(primary_hits, key=lambda x: x.q_st)
    
    # Get strand from best primary alignment (by alignment length)
    best_hit = max(hits, key=lambda x: x.blen)
    strand = best_hit.strand
    
    # Only use hits on the same strand
    hits = [h for h in hits if h.strand == strand]
    if not hits:
        return None
        
    # Find outermost bounds
    start = min(h.q_st for h in hits)
    end = max(h.q_en for h in hits)
    
    # Calculate coverage
    coverage = (end - start) / seq_len
    
    return start, end, coverage, strand

def process_sequence(aligner, name: str, seq: str, min_coverage=0.2):
    """
    Align sequence to reference, orient and clip if necessary.
    
    Args:
        aligner: Initialized minimap2 aligner
        name: Sequence name/header
        seq: Sequence content
        min_coverage: Minimum coverage threshold
        
    Returns:
        Tuple of (name, seq, description) or None if sequence doesn't meet criteria
    """
    # Get all alignments
    hits = list(aligner.map(seq))
    
    # Get alignment bounds and metrics
    alignment_info = get_alignment_bounds(hits)
    if not alignment_info:
        return None
        
    start, end, coverage, strand = alignment_info
    
    # Check coverage threshold
    if coverage < min_coverage:
        return None
    
    # Clip unaligned ends
    clipped_seq = seq[start:end]

    # Orient sequence if needed
    if strand == -1:
        clipped_seq = mp.revcomp(clipped_seq)
    
    # Create description with metrics
    description = f"coverage={coverage:.3f}"
    
    return name, clipped_seq, description

# Validate input files
input_fasta = snakemake.input['sequences']
ref_fasta = snakemake.input['reference']

if not os.path.exists(input_fasta):
    raise FileNotFoundError(f"Input sequences file {input_fasta} does not exist")
if not os.path.exists(ref_fasta):
    raise FileNotFoundError(f"Reference file {ref_fasta} does not exist")

output_fasta = snakemake.output[0]

# Get parameters with validation
min_coverage = snakemake.params.get("min_coverage", 0.2)
if not isinstance(min_coverage, (int, float)) or min_coverage <= 0 or min_coverage > 1:
    raise ValueError(f"min_coverage must be a number between 0 and 1, got {min_coverage}")

include_reference = snakemake.params.get("include_reference", True)
if not isinstance(include_reference, bool):
    raise ValueError(f"include_reference must be a boolean, got {include_reference}")

# Read reference sequence
try:
    ref_name, ref_seq, _, comment = next(mp.fastx_read(ref_fasta, read_comment=True))
except StopIteration:
    raise ValueError(f"Reference file {ref_fasta} is empty")
except Exception as e:
    raise ValueError(f"Failed to read reference file {ref_fasta}: {str(e)}")

# Initialize aligner
try:
    aligner = mp.Aligner(ref_fasta, preset="map-ont")
    if not aligner:
        raise Exception("Failed to initialize aligner")
except Exception as e:
    raise RuntimeError(f"Failed to initialize minimap2 aligner: {str(e)}")

# Process sequences
processed_records = []
try:
    for name, seq, _, comment in mp.fastx_read(input_fasta, read_comment=True):
        header = (name+' '+comment) if comment else name
        processed = process_sequence(
            aligner,
            header,
            seq,
            min_coverage=min_coverage
        )
        if processed:
            processed_records.append(processed)
except Exception as e:
    raise RuntimeError(f"Error processing sequences: {str(e)}")

# Write output
try:
    with open(output_fasta, 'w') as f:
        if include_reference:
            f.write(f">{ref_name}\n{ref_seq}\n")
        for name, seq, description in processed_records:
            f.write(f">{name} {description}\n{seq}\n")
except Exception as e:
    raise RuntimeError(f"Failed to write output file {output_fasta}: {str(e)}")

# Log summary
print(f"Processed {len(processed_records)} sequences")
print(f"Output written to {output_fasta}") 