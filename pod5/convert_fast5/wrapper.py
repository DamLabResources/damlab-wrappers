"""Wrapper for pod5 convert fast5.

This wrapper provides a Snakemake interface to the pod5 convert fast5 tool,
allowing users to convert Oxford Nanopore FAST5 files to the newer POD5 format.
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
fast5_files: str = ""

if 'fast5_files' in snakemake.input.keys():
    fast5_files = snakemake.input.fast5_files
elif 'fast5_dir' in snakemake.input.keys():
    fast5_dir = snakemake.input.fast5_dir
    # Validate input directory exists
    if not os.path.exists(fast5_dir):
        raise FileNotFoundError(f"Input directory {fast5_dir} does not exist")
    
    # Find files recursively in fast5_dir
    fast5_paths = list(Path(fast5_dir).rglob("*.fast5"))
    if not fast5_paths:
        raise ValueError(f"No FAST5 files found in directory {fast5_dir}")
    
    fast5_files = ' '.join(str(path) for path in fast5_paths)
else:
    raise ValueError("No input files or directories provided. Specify either 'fast5_files' or 'fast5_dir'")

# Validate output path
output_pod5 = snakemake.output[0]
output_dir = os.path.dirname(output_pod5)
if output_dir and not os.path.exists(output_dir):
    os.makedirs(output_dir, exist_ok=True)

# Get optional parameters
extra: str = snakemake.params.get("extra", "")
threads: Optional[int] = snakemake.threads
threads_arg: str = f"-t {threads}" if threads else ""

# Setup log
log: str = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Convert FAST5 to POD5
shell(
    f"pod5 convert fast5 {extra} --output {output_pod5} {threads_arg} {fast5_files} {log}"
) 