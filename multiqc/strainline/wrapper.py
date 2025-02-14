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

# Get input file
haplotypes_file = snakemake.input[0]
if not os.path.exists(haplotypes_file):
    raise ValueError("Input must be a haplotypes.final.fa file")

# Get output log file
output_log = snakemake.output[0]


def parse_haplotype_description(description):
    """Parse haplotype description string into dictionary."""
    #hap1 7973550x freq=1.000
    _, _, freq_str = description.split(' ')
    freq = float(freq_str.split('=')[1])
    return {
        'haplotype_freq': freq
    }

def get_haplotype_stats(haplotypes_file):
    """Get haplotype statistics from Strainline haplotypes file."""
    # Read haplotypes file
    with open(haplotypes_file, 'r') as f:
        for seqR in SeqIO.parse(f, 'fasta'):
            info = {}
            info['haplotype_id'] = seqR.id
            info.update(parse_haplotype_description(seqR.description))
            info['haplotype_length'] = len(seqR.seq)
            yield info


def parse_strainline_haplotypes(haplotypes_file, sample_name=None):
    """Parse Strainline haplotypes file for relevant metrics."""
    # Get sample name from parent directory name
    if sample_name is None:
        sample_name = os.path.basename(os.path.dirname(os.path.dirname(haplotypes_file)))

    haplotype_stats = pd.DataFrame(list(get_haplotype_stats(haplotypes_file)))

    metrics = {
        'sample_name': sample_name,
        'haplotype_count': len(haplotype_stats.index),
        'haplotype_max_length': int(haplotype_stats['haplotype_length'].max()),
        'haplotype_min_length': int(haplotype_stats['haplotype_length'].min()),
        'haplotype_mean_length': int(haplotype_stats['haplotype_length'].mean()),
        
        'haplotype_max_freq': float(haplotype_stats['haplotype_freq'].max()),
        'haplotype_min_freq': float(haplotype_stats['haplotype_freq'].min()),
        'haplotype_mean_freq': float(haplotype_stats['haplotype_freq'].mean()),
    }
    
    return metrics

def write_multiqc_log(metrics, output_file):
    """Write metrics in MultiQC-compatible YAML format."""
    with open(output_file, 'w') as handle:
        handle.write(f"# Strainline MultiQC Log\n")
        yaml.dump(metrics, handle, default_flow_style=False)

# Parse the haplotypes file
metrics = parse_strainline_haplotypes(haplotypes_file, sample_name=snakemake.params.get('sample_name'))

# Write the log file
write_multiqc_log(metrics, output_log) 