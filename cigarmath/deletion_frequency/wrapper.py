"""Wrapper for calculating deletion frequencies from BAM files"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2024"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"

import pysam # type: ignore
import yaml
import cigarmath as cm # type: ignore

if "snakemake" not in locals():
    import snakemake # type: ignore

# Get input/output files
input_bam = snakemake.input[0]
output_yaml = snakemake.output[0]

# Get required parameters
required_region = snakemake.params.get("required_region")
deletion_region = snakemake.params.get("deletion_region")
region_name = snakemake.params.get("region_name", "region")
sample_name = snakemake.params.get("sample_name", "sample")

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
    "deletion_frequency": 0.0
}

req_block_size = req_end - req_start
del_block_size = del_end - del_start
for start, cigars, segments in combined_segment_stream:
    if not cigars:
        continue
    stats["total_reads"] += 1
    ref_block = cm.reference_block(cigars, start)
    if cm.block_overlap_length(ref_block, (req_start, req_end)) == req_block_size:
        stats["reads_covering_required"] += 1
        for del_block in cm.reference_deletion_blocks(cigars, start, min_size=del_block_size    ):
            if cm.block_overlap_length(del_block, (del_start, del_end)) == del_block_size:
                stats["reads_with_deletion"] += 1
                break


try:
    print('did you make it here silly!')
    stats["deletion_frequency"] = stats["reads_with_deletion"] / stats["reads_covering_required"]
except ZeroDivisionError:
    stats["deletion_frequency"] = None


# Write output
with open(output_yaml, 'w') as f:
    f.write('# Cigarmath Deletion Frequency\n')
    yaml.dump(stats, f) 