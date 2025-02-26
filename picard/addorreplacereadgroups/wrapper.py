"""Wrapper for Picard AddOrReplaceReadGroups that accepts individual read group parameters"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2024"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"

import tempfile
from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts
import shlex

# Required read group parameters
REQUIRED_RG_PARAMS = ['ID', 'LB', 'PL', 'PU', 'SM']

# Optional read group parameters
OPTIONAL_RG_PARAMS = ['CN', 'DS', 'DT', 'FO', 'KS', 'PG', 'PI', 'PM']

if "snakemake" not in locals():
    #  Keeps linters happy but doesn't impact function
    import snakemake

# Get Java options
java_opts = get_java_opts(snakemake)

# Get logging
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Build read group arguments
rg_args = []

# Check required parameters
for param in REQUIRED_RG_PARAMS:
    value = snakemake.params.get(param)
    if value is None:
        raise ValueError(f"Required read group parameter '{param}' not provided")
    if param == 'ID':
        rg_args.append(f"--RGID {shlex.quote(str(value))}")
    else:
        rg_args.append(f"--RG{param} {shlex.quote(str(value))}")

# Add optional parameters if provided
for param in OPTIONAL_RG_PARAMS:
    value = snakemake.params.get(param)
    if value is not None:
        rg_args.append(f"--RG{param} {shlex.quote(str(value))}")

# Get any additional parameters
extra = snakemake.params.get("extra", "")

# Join all RG args
rg_args_str = " ".join(rg_args)

# Create command
with tempfile.TemporaryDirectory() as tmpdir:
    cmd = (
        "picard AddOrReplaceReadGroups"
        f" {java_opts}"
        f" --INPUT {snakemake.input[0]}"
        f" --OUTPUT {snakemake.output[0]}"
        f" --TMP_DIR {tmpdir}"
        f" {rg_args_str}"
        f" {extra}"
        f" {log}"
    )
    shell(cmd) 