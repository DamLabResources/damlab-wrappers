"""Wrapper for creating MultiQC-compatible logs from Strainline output"""

import os
from pathlib import Path

from snakemake.shell import shell # type: ignore
import yaml # type: ignore
from Bio import SeqIO # type: ignore
import pandas as pd # type: ignore

if "snakemake" not in locals():
    #  Keeps linters happy but doesn't impact funtion
    import snakemake # type: ignore

# Get input directory
input_dir = snakemake.input[0]
if not os.path.isdir(input_dir):
    raise ValueError("Input must be a Strainline output directory")

# Get output log file
output_log = snakemake.output[0]


def parse_haplotype_description(description):
    """Parse haplotype description string into dictionary."""
    # 7973550x freq=1.000
    _, freq_str = description.split(' ')
    freq = float(freq_str.split('=')[1])
    return {
        'haplotype_freq': freq
    }

def get_haplotype_stats(directory):
    """Get haplotype statistics from Strainline output directory."""
    haplotypes_file = os.path.join(directory, 'haplotypes.final.fa')
    if not os.path.exists(haplotypes_file):
        raise ValueError("Haplotypes file not found in Strainline output directory")
    
    # Read haplotypes file
    with open(haplotypes_file, 'r') as f:
        for seqR in SeqIO.parse(f, 'fasta'):
            info = {}
            info['haplotype_id'] = seqR.id
            info.update(parse_haplotype_description(seqR.description))
            info['haplotype_length'] = len(seqR.seq)
            yield info
    


def parse_strainline_dir(directory):
    """Parse Strainline output directory for relevant metrics."""

    haplotype_stats = pd.DataFrame(list(get_haplotype_stats(directory)))

    metrics = {
        'sample_id': os.path.basename(directory),
        'haplotype_count': len(haplotype_stats.index),
        'haplotype_max_length': haplotype_stats['haplotype_length'].max(),
        'haplotype_min_length': haplotype_stats['haplotype_length'].min(),
        'haplotype_mean_length': haplotype_stats['haplotype_length'].mean(),
        
        'haplotype_max_freq': haplotype_stats['haplotype_freq'].max(),
        'haplotype_min_freq': haplotype_stats['haplotype_freq'].min(),
        'haplotype_mean_freq': haplotype_stats['haplotype_freq'].mean(),
    }
    
    
    return metrics

def write_multiqc_log(metrics, output_file):
    """Write metrics in MultiQC-compatible YAML format."""
    with open(output_file, 'w') as handle:
        handle.write(f"# Strainline MultiQC Log\n")
        yaml.dump(metrics, handle, default_flow_style=False)

# Parse the directory
metrics = parse_strainline_dir(input_dir)

# Write the log file
write_multiqc_log(metrics, output_log) 