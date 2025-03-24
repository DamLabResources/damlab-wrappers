"""Wrapper for calculating HIV sequence intactness statistics"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2024"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"

import pysam
from Bio import SeqIO
import yaml
from pathlib import Path

if "snakemake" not in locals():
    import snakemake # type: ignore


def stream_reads(path, chrom):

    # Determine file type from extension
    suffix = Path(path).suffix.lower()

    if suffix in ['.bam', '.sam']:
        # Process BAM/SAM file
        with pysam.AlignmentFile(path, 'r') as f:
            for read in f:
                if chrom is not None:
                    if read.reference_name != chrom:
                        continue
                if read.is_unmapped:
                    seq_len = len(read.query_sequence)
                else:
                    seq_len = read.query_alignment_length
                yield seq_len

    elif suffix in ['.fastq', '.fasta']:
        # Process FASTA/FASTQ file
        with open(path, 'r') as f:
            for record in SeqIO.parse(f, 'fastq' if suffix == '.fastq' else 'fasta'):
                yield len(record.seq)


def count_sequences_by_length(path, min_countable, max_countable, min_intact, max_intact, chrom):
    """Count sequences within countable and intact length ranges"""
    stats = {
        'total_reads': 0,
        'countable_reads': 0,
        'intact_reads': 0
    }
    
    # Determine file type from extension
    suffix = Path(path).suffix.lower()

    for seq_len in stream_reads(path, chrom):
        stats['total_reads'] += 1
        if min_countable <= seq_len <= max_countable:
            stats['countable_reads'] += 1
        if min_intact <= seq_len <= max_intact:
            stats['intact_reads'] += 1
    
    # Calculate percentages
    try:
        stats['percent_countable'] = (stats['countable_reads'] / stats['total_reads']) * 100
    except ZeroDivisionError:
        stats['percent_countable'] = 0

    try:
        stats['percent_intact'] = (stats['intact_reads'] / stats['countable_reads']) * 100
    except ZeroDivisionError:
        stats['percent_intact'] = 0
        
    return stats

# Get input/output paths
input_path = str(snakemake.input[0])
output_path = str(snakemake.output[0])

# Get parameters
min_countable = snakemake.params['min_countable']
max_countable = snakemake.params['max_countable']
min_intact = snakemake.params['min_intact']
max_intact = snakemake.params['max_intact']
sample_name = snakemake.params.get('sample_name', 'sample')
chrom = snakemake.params.get('chrom')

# Calculate statistics
stats = count_sequences_by_length(
    input_path,
    min_countable,
    max_countable,
    min_intact,
    max_intact,
    chrom
)

# Add sample name if provided
stats['sample_name'] = sample_name

# Write output
with open(output_path, 'w') as f:
    f.write('# HIV Sequence Intactness Statistics\n')
    yaml.dump(stats, f) 