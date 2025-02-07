"""Wrapper for pod5 subset by channel"""

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
if 'pod5_files' in snakemake.input.keys():
    pod5_files = snakemake.input.pod5_files
elif 'pod5_dir' in snakemake.input.keys():
    pod5_dir = snakemake.input.pod5_dir
    # Find files recursively in pod5_dir
    pod5_files = ' '.join(str(path) for path in Path(pod5_dir).rglob("*.pod5"))
    
else:
    raise ValueError("No input files or directories provided")

output_dir = snakemake.output[0]
Path(output_dir).mkdir(parents=True, exist_ok=True)

# Get optional parameters
extra = snakemake.params.get("extra", "")

read2channel_file = Path(output_dir) / "read2channel.tsv"

# Setup log
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

read2channel_cmd = f"pod5 view -r {pod5_files} --include 'read_id,channel' --output {read2channel_file} {log}"

# Create read2channel file
shell(read2channel_cmd)

subset_cmd = f"pod5 subset -r {pod5_files} --summary {read2channel_file} --columns channel --output {output_dir} {extra} {log}"

# Split by channel
shell(subset_cmd)
