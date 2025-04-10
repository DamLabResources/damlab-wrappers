"""Wrapper for pod5 subset by channel.

This wrapper provides a Snakemake interface to the pod5 subset tool,
allowing users to split POD5 files by channel for downstream processing.
"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2024"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.0.0"

import os
from pathlib import Path
from typing import List, Optional, Union

from snakemake.shell import shell  # type: ignore

# This is a common pattern in Snakemake wrappers
# It allows the wrapper to be imported without snakemake being in the global namespace
# This is useful for testing and linting
if "snakemake" not in locals():
    import snakemake  # type: ignore

# Check if version is specified and compatible
if hasattr(snakemake, "params") and "version" in snakemake.params:
    requested_version = snakemake.params.version
    if requested_version != __version__:
        print(f"Warning: Requested version {requested_version} does not match wrapper version {__version__}")

# Extract arguments from snakemake object
pod5_files: str = ""

if 'pod5_files' in snakemake.input.keys():
    pod5_files = snakemake.input.pod5_files
elif 'pod5_dir' in snakemake.input.keys():
    pod5_dir = snakemake.input.pod5_dir
    # Validate input directory exists
    if not os.path.exists(pod5_dir):
        raise FileNotFoundError(f"Input directory {pod5_dir} does not exist")
    
    # Find files recursively in pod5_dir
    pod5_paths = list(Path(pod5_dir).rglob("*.pod5"))
    if not pod5_paths:
        raise ValueError(f"No POD5 files found in directory {pod5_dir}")
    
    pod5_files = ' '.join(str(path) for path in pod5_paths)
else:
    raise ValueError("No input files or directories provided. Specify either 'pod5_files' or 'pod5_dir'")

# Validate output directory
output_dir = snakemake.output[0]
Path(output_dir).mkdir(parents=True, exist_ok=True)

# Get optional parameters
extra: str = snakemake.params.get("extra", "")

# Create read2channel file path
read2channel_file = Path(output_dir) / "read2channel.tsv"

# Setup log
log: str = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Create read2channel file
read2channel_cmd = f"pod5 view -r {pod5_files} --include 'read_id,channel' --output {read2channel_file} {log}"
shell(read2channel_cmd)

# Split by channel
subset_cmd = f"pod5 subset -r {pod5_files} --summary {read2channel_file} --columns channel --output {output_dir} {extra} {log}"
shell(subset_cmd)
