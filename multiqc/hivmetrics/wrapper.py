"""Wrapper for creating MultiQC-compatible logs from HIV sequencing data"""

import os
from pathlib import Path
import pysam
import numpy as np
import yaml
from snakemake.shell import shell

if "snakemake" not in locals():
    import snakemake # type: ignore

def create_length_histogram(lengths, max_size=10000, bin_size=100):
    """Create histogram of read lengths."""
    bins = np.arange(0, max_size + bin_size, bin_size)
    hist, edges = np.histogram(lengths, bins=bins)
    return {
        'counts': hist.tolist(),
        'bin_edges': edges.tolist(),
        'bin_size': bin_size
    }

def parse_hiv_bam(bam_file, max_size=10000, sample_name=None):
    """Parse BAM file for HIV metrics."""
    if sample_name is None:
        sample_name = Path(bam_file).stem

    metrics = {
        'sample_name': sample_name,
        'total_reads': 0,
        'mapped_reads': 0,
        'supplementary_reads': 0,
        'mapped_lengths': [],
        'length_histogram': None
    }
    
    # Read BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            metrics['total_reads'] += 1
            
            if read.is_supplementary:
                metrics['supplementary_reads'] += 1
            
            if not read.is_unmapped:
                metrics['mapped_reads'] += 1
                metrics['mapped_lengths'].append(read.query_length)
    
    # Create length histogram
    if metrics['mapped_lengths']:
        metrics['length_histogram'] = create_length_histogram(
            metrics['mapped_lengths'], 
            max_size=max_size
        )
        
        # Add summary statistics
        metrics['mean_mapped_length'] = float(np.mean(metrics['mapped_lengths']))
        metrics['median_mapped_length'] = float(np.median(metrics['mapped_lengths']))
        
    else:
        metrics['mean_mapped_length'] = 0
        metrics['median_mapped_length'] = 0
    
    # Calculate mapping rate
    metrics['mapping_rate'] = metrics['mapped_reads'] / metrics['total_reads'] if metrics['total_reads'] > 0 else 0
    
    # Clean up temporary data
    del metrics['mapped_lengths']
    
    return metrics

def write_multiqc_log(metrics, output_file):
    """Write metrics in MultiQC-compatible YAML format."""
    with open(output_file, 'w') as handle:
        handle.write(f"# HIVmetrics MultiQC Log\n")
        yaml.dump(metrics, handle, default_flow_style=False)

# Get input/output files
bam_file = snakemake.input[0]
output_log = snakemake.output[0]

if not os.path.exists(bam_file):
    raise ValueError("Input must be a BAM file")

# Get parameters
max_size = snakemake.params.get('max_size', 10000)

# Parse the BAM file
metrics = parse_hiv_bam(
    bam_file,
    max_size=max_size,
    sample_name=snakemake.params.get('sample_name')
)

# Write the log file
write_multiqc_log(metrics, output_log) 