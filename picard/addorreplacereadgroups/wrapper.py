"""Wrapper for Picard AddOrReplaceReadGroups.

This wrapper provides a Snakemake interface to Picard's AddOrReplaceReadGroups tool,
allowing users to add or replace read groups in BAM files with a structured parameter interface.
"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2024"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.0.0"

import os
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Union

from snakemake.shell import shell  # type: ignore
from snakemake_wrapper_utils.java import get_java_opts

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

# Required read group parameters
REQUIRED_RG_PARAMS: List[str] = ['ID', 'LB', 'PL', 'PU', 'SM']

# Optional read group parameters
OPTIONAL_RG_PARAMS: List[str] = ['CN', 'DS', 'DT', 'FO', 'KS', 'PG', 'PI', 'PM']

# Get input and output files
input_bam: str = snakemake.input[0]
output_bam: str = snakemake.output[0]

# Validate input file exists
if not os.path.exists(input_bam):
    raise FileNotFoundError(f"Input file {input_bam} does not exist")

# Get Java options
java_opts: str = get_java_opts(snakemake)

# Get logging
log: str = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Build read group arguments
rg_args: List[str] = []

# Check required parameters
for param in REQUIRED_RG_PARAMS:
    value = snakemake.params.get(param)
    if value is None:
        raise ValueError(f"Required read group parameter '{param}' not provided")
    if param == 'ID':
        rg_args.append(f"--RGID {value}")
    else:
        rg_args.append(f"--RG{param} {value}")

# Add optional parameters if provided
for param in OPTIONAL_RG_PARAMS:
    value = snakemake.params.get(param)
    if value is not None:
        rg_args.append(f"--RG{param} {value}")

# Get any additional parameters
extra: str = snakemake.params.get("extra", "")
validation_stringency: Optional[str] = snakemake.params.get("validation_stringency", None)
create_index: bool = snakemake.params.get("create_index", True)
compression_level: Optional[int] = snakemake.params.get("compression_level", None)

# Add validation stringency if specified
if validation_stringency:
    rg_args.append(f"--VALIDATION_STRINGENCY {validation_stringency}")

# Add create index parameter
rg_args.append(f"--CREATE_INDEX {str(create_index).lower()}")

# Add compression level if specified
if compression_level is not None:
    rg_args.append(f"--COMPRESSION_LEVEL {compression_level}")

# Join all RG args
rg_args_str: str = " ".join(rg_args)

# Create command
with tempfile.TemporaryDirectory() as tmpdir:
    cmd = (
        "picard AddOrReplaceReadGroups"
        f" {java_opts}"
        f" --INPUT {input_bam}"
        f" --OUTPUT {output_bam}"
        f" --TMP_DIR {tmpdir}"
        f" {rg_args_str}"
        f" {extra}"
        f" {log}"
    )
    shell(cmd) 