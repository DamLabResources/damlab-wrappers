"""Wrapper for seqkit amplicon to check primers"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2024"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"

from snakemake.shell import shell # type: ignore
import csv
from pathlib import Path
from collections import defaultdict
import yaml
from itertools import combinations

if "snakemake" not in locals():
    #  Keeps linters happy but doesn't impact function
    import snakemake # type: ignore

# Extract arguments from snakemake object
reads = snakemake.input.reads
primers = snakemake.input.primers
output_csv = snakemake.output[0]

# Get optional summary output if specified
summary_yaml = snakemake.output.get("summary", None)

# Get optional parameters
extra = snakemake.params.get("extra", "")

# Create temporary directory for intermediate files
temp_dir = Path("temp_primercheck")
temp_dir.mkdir(exist_ok=True)
temp_bed = temp_dir / "amplicons.bed"

# Setup log
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Run seqkit amplicon
shell(
    f"seqkit amplicon --bed --primer-file {primers} {reads} > {temp_bed} {log}"
)

# Read BED file and create summary
# BED format: chrom start end name score strand amplicon_length
results = defaultdict(dict)
primer_names = set()

with open(temp_bed) as f:
    for line in f:
        fields = line.strip().split('\t')
        read_id = fields[0]
        primer_name = fields[3]
        amplicon_length = len(fields[6])  # Length of the amplicon sequence
        
        results[read_id][primer_name] = amplicon_length
        primer_names.add(primer_name)

# Convert to CSV
with open(output_csv, 'w', newline='') as f:
    writer = csv.writer(f)
    
    # Write header
    header = ['read_id'] + sorted(list(primer_names))
    writer.writerow(header)
    
    # Write data
    for read_id in sorted(results.keys()):
        row = [read_id]
        for primer in header[1:]:  # Skip read_id column
            row.append(results[read_id].get(primer, 'None'))
        writer.writerow(row)

# Generate summary statistics if requested
if summary_yaml:
    # Count total sequences checked
    total_seqs = len(results)
    
    # Count hits per primer
    primer_hits = {primer: sum(1 for read in results.values() if primer in read) 
                  for primer in primer_names}
    
    # Generate pairwise matrix
    pairwise_hits = {}
    for p1, p2 in combinations(sorted(primer_names), 2):
        dual_hits = sum(1 for read in results.values() 
                       if p1 in read and p2 in read)
        pairwise_hits[f"{p1}_x_{p2}"] = dual_hits

    summary = {
        'total_sequences': total_seqs,
        'primer_hits': primer_hits,
        'pairwise_hits': pairwise_hits
    }
    
    with open(summary_yaml, 'w') as f:
        f.write("# Seqkit Primer Check Summary\n")
        yaml.dump(summary, f, default_flow_style=False)

# Cleanup
shell(f"rm -rf {temp_dir}") 