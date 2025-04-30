"""Wrapper for creating gel visualizations from BAM files using cigarmath"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2025"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.0.0"

import yaml
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import cigarmath as cm
from typing import Dict, Any, Generator
import numpy as np
import logging

# Configure logging
if "snakemake" in locals():
    log_file = snakemake.log[0] if hasattr(snakemake, 'log') and snakemake.log else None
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        filename=log_file
    )
else:
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

logger = logging.getLogger('bam2gel-wrapper')

if "snakemake" not in locals():
    import snakemake

# Get input/output files
input_bams = snakemake.input
output_png = snakemake.output[0]
output_yaml = snakemake.output[1] if len(snakemake.output) > 1 else None

# Get parameters
params = snakemake.params
mode = params.get("mode", "read_length")
names = params.get("names", [f"Sample_{i+1}" for i in range(len(input_bams))])
bin_width = params.get("bin_width", 10)
max_size = params.get("max_size", 1000)

# Get plot parameters
plot_kwargs = params.get("plot_kwargs", {})
fig_kwargs = params.get("fig_kwargs", {"figsize": (10, 6)})
axis_kwargs = params.get("axis_kwargs", {})

def process_bam(bam_file: str, name: str) -> Generator[Dict[str, Any], None, None]:
    """Process a single BAM file and yield size information based on mode"""
    
    segment_stream = cm.io.segment_stream_pysam(bam_file, mode='rb')
    combined_segment_stream = cm.io.combined_segment_stream(segment_stream)
    
    for start, cigars, segments in combined_segment_stream:
        if not cigars:
            continue
            
        if mode == "read_length":
            size = segments[0].query_length
        elif mode == "reference_block_size":
            ref_block = cm.reference_block(cigars, start)
            size = ref_block[1] - ref_block[0]
        elif mode == "query_block_size":
            query_block = cm.query_block(cigars)
            size = query_block[1] - query_block[0]
        else:
            raise ValueError(f"Invalid mode: {mode}")
            
        if size <= max_size:
            yield {"size": size, "lane": name}

# Process all BAM files and collect data
all_data = []
metrics = {"mode": mode, "lanes": {}}

for bam_file, name in zip(input_bams, names):
    lane_data = []
    for item in process_bam(bam_file, name):
        lane_data.append(item)
    
    # Calculate lane metrics and convert numpy types to native Python types
    lane_metrics = {
        "total_reads": int(len(lane_data)),
        "mean_size": float(np.mean([d["size"] for d in lane_data])) if lane_data else 0.0,
        "median_size": float(np.median([d["size"] for d in lane_data])) if lane_data else 0.0,
        "min_size": int(min([d["size"] for d in lane_data])) if lane_data else 0,
        "max_size": int(max([d["size"] for d in lane_data])) if lane_data else 0
    }
    metrics["lanes"][name] = lane_metrics
    all_data.extend(lane_data)

# Create DataFrame
df = pd.DataFrame(all_data)

# Create gel visualization
plt.figure(**fig_kwargs)

# Set default plot kwargs if not provided
default_plot_kwargs = {
    "x": "lane",
    "y": "size",
    "bins": (len(names), np.arange(0, max_size + bin_width, bin_width)),
    "palette": "husl"
}
plot_kwargs = {**default_plot_kwargs, **plot_kwargs}

# Create the plot
sns.histplot(data=df, **plot_kwargs)

# Set default axis kwargs if not provided
default_axis_kwargs = {
    "xlabel": "Lane",
    "ylabel": "Size (bp)",
    "title": f"Gel Visualization - {mode.replace('_', ' ').title()}"
}
axis_kwargs = {**default_axis_kwargs, **axis_kwargs}

# Apply axis customizations
for key, value in axis_kwargs.items():
    if hasattr(plt, key):
        getattr(plt, key)(value)
    elif hasattr(plt, f"set_{key}"):
        getattr(plt, f"set_{key}")(value)

# Save plot
plt.savefig(output_png, dpi=300, bbox_inches='tight')
plt.close()

# Save metrics if requested
if output_yaml:
    with open(output_yaml, 'w') as f:
        f.write('# Cigarmath BAM2GEL Metrics\n')
        yaml.dump(metrics, f) 