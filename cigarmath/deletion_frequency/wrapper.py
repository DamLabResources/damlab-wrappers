"""Wrapper for calculating deletion frequencies from BAM files"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2024"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.0.0"

import pysam # type: ignore
import yaml
import csv
import cigarmath as cm # type: ignore

if "snakemake" not in locals():
    import snakemake # type: ignore

# Get input/output files
input_bam = snakemake.input[0]
output_yaml = snakemake.output[0]
# Check if CSV output is specified
output_csv = snakemake.output[1] if len(snakemake.output) > 1 else None

# Get required parameters
required_region = snakemake.params.get("required_region")
deletion_region = snakemake.params.get("deletion_region")
region_name = snakemake.params.get("region_name", "region")
sample_name = snakemake.params.get("sample_name", "sample")
min_deletion_size = snakemake.params.get("min_deletion_size", 1)

# Validate parameters
if not required_region or not deletion_region:
    raise ValueError("Both required_region and deletion_region parameters must be specified")

# Parse region strings (format: "chr:start-end")
try:
    req_chrom, req_coords = required_region.split(":")
    req_start, req_end = map(int, req_coords.split("-"))
    
    del_chrom, del_coords = deletion_region.split(":")
    del_start, del_end = map(int, del_coords.split("-"))
except ValueError:
    raise ValueError("Region format must be 'chr:start-end'")

segment_stream = cm.io.segment_stream_pysam(input_bam, mode='rb')
fixed_segment_stream = (segment for segment in segment_stream if segment.cigartuples)
combined_segment_stream = cm.io.combined_segment_stream(fixed_segment_stream)

stats = {
    "sample_name": sample_name,
    "region_name": region_name,
    "total_reads": 0,
    "reads_covering_required": 0,
    "reads_with_deletion": 0,
    "reads_covered_with_deletion": 0,
    "deletion_frequency": 0.0
}

# List to store read-level information if CSV output is requested
read_data = []

req_block_size = req_end - req_start
del_block_size = del_end - del_start
for start, cigars, segments in combined_segment_stream:
    if not cigars:
        continue
    
    stats["total_reads"] += 1
    
    # Get read name from the first segment
    read_name = segments[0].query_name if segments else "unknown"
    
    # Check if read covers required region
    ref_block = cm.reference_block(cigars, start)
    covers_required = cm.block_overlap_length(ref_block, (req_start, req_end)) == req_block_size
    
    if covers_required:
        stats["reads_covering_required"] += 1
    
    # Check if read has deletion in specified region
    has_deletion = False
    for del_block in cm.reference_deletion_blocks(cigars, start, min_size=min_deletion_size):
        if del_block[1] < del_start:
            # Deletion is before the start of the region, so we can skip
            continue
        elif del_block[0] > del_end:
            # All deletions are past the end of the region, so we can break
            break
        elif cm.block_overlap_length(del_block, (del_start, del_end)) > min_deletion_size:
            has_deletion = True
            break
    
    # Only increment the counter if the read has a deletion
    if has_deletion:
        stats["reads_with_deletion"] += 1

    stats["reads_covered_with_deletion"] += has_deletion & covers_required
    
    # Store read-level information if CSV output is requested
    if output_csv:
        read_data.append({
            "read_name": read_name,
            "covers_required": covers_required,
            "has_deletion": has_deletion
        })

try:
    stats["deletion_frequency"] = stats["reads_covered_with_deletion"] / stats["reads_covering_required"]
except ZeroDivisionError:
    stats["deletion_frequency"] = None

# Write YAML output
with open(output_yaml, 'w') as f:
    f.write('# Cigarmath Deletion Frequency\n')
    yaml.dump(stats, f)

# Write CSV output if requested
if output_csv:
    with open(output_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["read_name", "covers_required", "has_deletion"])
        writer.writeheader()
        writer.writerows(read_data) 