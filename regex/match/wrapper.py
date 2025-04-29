"""Wrapper for regex pattern matching in FASTA files"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2025"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.0.0"

import regex
import pandas as pd
from Bio import SeqIO
import yaml


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

def process_fasta(fasta_path, patterns):
    """Process FASTA file and extract matches for each pattern"""
    results = []
    
    for record in SeqIO.parse(fasta_path, "fasta"):
        row = {'read_name': record.id}
        
        # Extract matches for each pattern
        for i, pattern in enumerate(patterns):
            match = extract_matches(pattern, str(record.seq))
            row[f'pattern_{i+1}'] = match
        
        results.append(row)
    
    return pd.DataFrame(results)

# Get input/output paths
fasta_in_path = str(snakemake.input[0])
csv_out_path = str(snakemake.output[0])
metrics_file = snakemake.output.get('metrics', None)

# Get regex patterns
patterns = get_regex_patterns(dict(snakemake.params))

# Process FASTA file
df = process_fasta(fasta_in_path, patterns)

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