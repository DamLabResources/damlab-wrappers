"""Wrapper for calculating pileup depth from BAM files"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2024"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"

import pysam # type: ignore
import cigarmath as cm # type: ignore
import csv
import math
from collections import defaultdict, Counter

if "snakemake" not in locals():
    import snakemake # type: ignore

def calculate_entropy(base_counts):
    """Calculate Shannon entropy of base frequencies at a position.
    
    Args:
        base_counts: Counter object with base counts
        
    Returns:
        float: Shannon entropy value
    """
    total = sum(base_counts.values())
    if total == 0:
        return 0.0
        
    entropy = 0.0
    for count in base_counts.values():
        if count == 0:
            continue
        freq = count / total
        entropy -= freq * math.log2(freq)
    return entropy

# Get input/output files
input_bam = snakemake.input[0]
output_tsv = snakemake.output[0]

# Get optional parameters
region = snakemake.params.get("region", None)
min_mapq = snakemake.params.get("min_mapq", 0)

# Initialize combined pileup counts
combined_counts = defaultdict(Counter)

# Stream segments from BAM file
segment_stream = cm.io.segment_stream_pysam(
    input_bam, 
    mode='rb',
    fetch=region,
    min_mapq=min_mapq
)

# Process each segment
for segment in segment_stream:
    if not segment.cigartuples or not segment.query_sequence:
        continue
        
    # Calculate depth for this segment
    counts = cm.depth(
        segment.cigartuples,
        reference_start=segment.reference_start,
        query_sequence=segment.query_sequence
    )
    
    # Update combined counts
    for pos, base_counts in counts.items():
        combined_counts[pos].update(base_counts)

# Define column order for bases (including gap)
BASES = ['A', 'C', 'G', 'T', 'N', '-']

# Write output TSV
with open(output_tsv, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    
    # Write header
    writer.writerow(['Position', 'Depth', 'Entropy'] + BASES)
    
    # Write counts for each position
    for pos in sorted(combined_counts.keys()):
        base_counts = combined_counts[pos]
        
        # Calculate non-gap depth
        depth = sum(count for base, count in base_counts.items() if base != '-')
        
        # Calculate entropy
        entropy = calculate_entropy(base_counts)
        
        # Get counts for each base (0 if not present)
        base_values = [base_counts[base] for base in BASES]
        
        # Write line
        writer.writerow([pos, depth, f'{entropy:.3f}'] + base_values) 