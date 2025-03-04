"""Wrapper for creating MultiQC-compatible logs from Dorado output"""

import os
from pathlib import Path
import pysam
import numpy as np
import yaml
from snakemake.shell import shell

if "snakemake" not in locals():
    import snakemake

def calculate_n50(lengths):
    """Calculate N50 from list of lengths."""
    sorted_lengths = sorted(lengths, reverse=True)
    total = sum(sorted_lengths)
    running_sum = 0
    for length in sorted_lengths:
        running_sum += length
        if running_sum >= total/2:
            return length
    return 0

def parse_dorado_bam(bam_file, sample_name=None):
    """Parse Dorado BAM file for relevant metrics."""
    if sample_name is None:
        sample_name = Path(bam_file).stem

    # Initialize storage for metrics
    metrics = {
        'sample_name': sample_name,
        'simplex': {
            'total_reads': 0,
            'total_bases': 0,
            'passed_reads': 0,
            'failed_reads': 0,
            'lengths': [],
            'qscores': [],
            'no_duplex': 0,  # dx:i:0
            'has_duplex': 0  # dx:i:-1
        },
        'duplex': {
            'total_reads': 0,
            'total_bases': 0,
            'passed_reads': 0,
            'failed_reads': 0,
            'lengths': [],
            'qscores': []
        }
    }
    
    # Read BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            # Get read type from dx tag
            dx_tag = read.get_tag('dx') if read.has_tag('dx') else 0
            
            # Determine read category
            if dx_tag == 1:
                category = 'duplex'
            else:
                category = 'simplex'
                if dx_tag == 0:
                    metrics['simplex']['no_duplex'] += 1
                elif dx_tag == -1:
                    metrics['simplex']['has_duplex'] += 1
            
            # Update metrics for this category
            metrics[category]['total_reads'] += 1
            metrics[category]['total_bases'] += read.query_length
            
            if read.is_qcfail:
                metrics[category]['failed_reads'] += 1
            else:
                metrics[category]['passed_reads'] += 1
            
            metrics[category]['lengths'].append(read.query_length)
            metrics[category]['qscores'].append(read.mapping_quality if read.mapping_quality is not None else 0)
    
    # Calculate summary statistics
    for category in ['simplex', 'duplex']:
        if metrics[category]['total_reads'] > 0:
            metrics[category].update({
                'mean_read_length': float(np.mean(metrics[category]['lengths'])),
                'read_length_n50': calculate_n50(metrics[category]['lengths']),
                'mean_qscore': float(np.mean(metrics[category]['qscores'])),
                'pass_rate': metrics[category]['passed_reads'] / metrics[category]['total_reads']
            })
        else:
            metrics[category].update({
                'mean_read_length': 0,
                'read_length_n50': 0,
                'mean_qscore': 0,
                'pass_rate': 0
            })
        
        # Clean up temporary lists
        del metrics[category]['lengths']
        del metrics[category]['qscores']
    
    return metrics

def write_multiqc_log(metrics, output_file):
    """Write metrics in MultiQC-compatible YAML format."""
    with open(output_file, 'w') as handle:
        handle.write(f"# Dorado MultiQC Log\n")
        yaml.dump(metrics, handle, default_flow_style=False)

# Get input/output files
bam_file = snakemake.input[0]
output_log = snakemake.output[0]

if not os.path.exists(bam_file):
    raise ValueError("Input must be a BAM file")

# Parse the BAM file
metrics = parse_dorado_bam(
    bam_file, 
    sample_name=snakemake.params.get('sample_name')
)

# Write the log file
write_multiqc_log(metrics, output_log) 