"""Wrapper for pod5 convert fast5"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2024"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"

from snakemake.shell import shell # type: ignore
from pathlib import Path

if "snakemake" not in locals():
    #  Keeps linters happy but doesn't impact function
    import snakemake # type: ignore

# Extract arguments from snakemake object

if 'fast5_files' in snakemake.input.keys():
    fast5_files = snakemake.input.fast5_files
elif 'fast5_dir' in snakemake.input.keys():
    fast5_dir = snakemake.input.fast5_dir   
    # Find files recursively in fast5_dir
    fast5_files = ' '.join(str(path) for path in Path(fast5_dir).rglob("*.fast5"))
else:
    raise ValueError("No input files or directories provided")


output_pod5 = snakemake.output[0]

# Get optional parameters
extra = snakemake.params.get("extra", "")

threads = f"-t {snakemake.threads}" if snakemake.threads else ""

# Setup log
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Convert FAST5 to POD5
shell(
    f"pod5 convert fast5 {extra} --output {output_pod5} {threads} {fast5_files} {log}"
) 