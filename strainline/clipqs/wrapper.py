"""Wrapper for clipping and orienting sequences using minimap2"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2024"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"

from typing import List, Optional, Tuple
import mappy as mp # type: ignore

if "snakemake" not in locals():
    import snakemake  # type: ignore

def get_alignment_bounds(hits: List[mp.Alignment]) -> Optional[Tuple[int, int, float, int]]:
    """
    Get the outermost alignment bounds from all hits
    Returns (start, end, coverage, strand) or None if no valid alignments
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
    Align sequence to reference, orient and clip if necessary
    Returns tuple of (name, seq, description) or None if sequence doesn't meet criteria
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
    
    # Orient sequence if needed
    if strand == -1:
        seq = mp.revcomp(seq)
    
    # Clip unaligned ends
    clipped_seq = seq[start:end]
    
    # Create description with metrics
    description = f"coverage={coverage:.3f}"
    
    return name, clipped_seq, description

# Get input/output files
input_fasta = snakemake.input['sequences']
ref_fasta = snakemake.input['reference']
output_fasta = snakemake.output[0]

# Get parameters
min_coverage = snakemake.params.get("min_coverage", 0.2)
include_reference = snakemake.params.get("include_reference", True)

ref_name, ref_seq, _, comment= next(mp.fastx_read(ref_fasta, read_comment=True))

# Initialize aligner
aligner = mp.Aligner(ref_fasta, preset="map-ont")
if not aligner:
    raise Exception("Failed to load reference")

# Process sequences
processed_records = []
for name, seq, _, comment in mp.fastx_read(input_fasta, read_comment=True):
    processed = process_sequence(
        aligner,
        name+' '+comment,
        seq,
        min_coverage=min_coverage
    )
    if processed:
        processed_records.append(processed)

# Write output
with open(output_fasta, 'w') as f:
    if include_reference:
        f.write(f">{ref_name}\n{ref_seq}\n")
    for name, seq, description in processed_records:
        f.write(f">{name} {description}\n{seq}\n") 