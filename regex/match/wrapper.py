"""Wrapper for regex pattern matching in FASTA/FASTQ files"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2024"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.0.0"

import regex
import pandas as pd
from Bio import SeqIO
import yaml
import gzip
import os

if "snakemake" not in locals():
    import snakemake # type: ignore

def get_regex_patterns(params):
    """Get regex patterns from parameters"""
    if 'patterns' not in params:
        raise ValueError("No regex patterns provided in params. Please provide patterns as a list in the 'patterns' parameter.")
    
    patterns = params['patterns']
    if not isinstance(patterns, list):
        raise ValueError("Patterns must be provided as a list")
    
    return [regex.compile(pattern, flags=regex.BESTMATCH+regex.IGNORECASE) for pattern in patterns]

def extract_matches(reg, text):
    """Extract first match from text using regex pattern"""
    match = reg.findall(text)
    if match:
        return match[0]
    return None

def get_sequence_reader(file_path):
    """Determine the appropriate sequence reader based on file extension"""
    _, ext = os.path.splitext(file_path)
    ext = ext.lower()
    
    # Handle gzipped files
    if ext == '.gz':
        _, base_ext = os.path.splitext(file_path[:-3])  # Remove .gz and get base extension
        base_ext = base_ext.lower()
        if base_ext in ['.fasta', '.fa', '.fna']:
            return lambda f: SeqIO.parse(gzip.open(f, 'rt'), 'fasta')
        elif base_ext in ['.fastq', '.fq']:
            return lambda f: SeqIO.parse(gzip.open(f, 'rt'), 'fastq')
        else:
            raise ValueError(f"Unsupported gzipped file format: {base_ext}")
    
    # Handle uncompressed files
    if ext in ['.fasta', '.fa', '.fna']:
        return lambda f: SeqIO.parse(f, 'fasta')
    elif ext in ['.fastq', '.fq']:
        return lambda f: SeqIO.parse(f, 'fastq')
    else:
        raise ValueError(f"Unsupported file format: {ext}")

def process_sequences(file_path, patterns):
    """Process sequence file and extract matches for each pattern"""
    results = []
    reader = get_sequence_reader(file_path)
    
    with open(file_path, 'rb') if not file_path.endswith('.gz') else gzip.open(file_path, 'rt') as handle:
        for record in reader(file_path):
            row = {'read_name': record.id}
            
            # Extract matches for each pattern
            for i, pattern in enumerate(patterns):
                match = extract_matches(pattern, str(record.seq))
                row[f'pattern_{i+1}'] = match
            
            results.append(row)
    
    return pd.DataFrame(results)

# Get input/output paths
seq_in_path = str(snakemake.input[0])
csv_out_path = str(snakemake.output[0])
metrics_file = snakemake.output.get('metrics', None)

# Get regex patterns
patterns = get_regex_patterns(dict(snakemake.params))

# Process sequence file
df = process_sequences(seq_in_path, patterns)

# Write results to CSV
df.to_csv(csv_out_path, index=False)

# Generate and write metrics if requested
if metrics_file:
    metrics = {
        'total_reads': int(len(df)),  # Convert to native Python int
        'pattern_matches': {
            f'pattern_{i+1}': int(df[f'pattern_{i+1}'].notna().sum())  # Convert to native Python int
            for i in range(len(patterns))
        }
    }
    
    with open(metrics_file, 'w') as f:
        f.write('# Regex pattern matching metrics\n')
        yaml.dump(metrics, f, default_flow_style=False) 