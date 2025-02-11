"""Wrapper for MultiQC with DAMlab plugins"""

from snakemake.shell import shell # type: ignore
import os

if "snakemake" not in locals():
    #  Keeps linters happy but doesn't impact funtion
    import snakemake # type: ignore

# Get input and output files
input_files = snakemake.input
output_report = snakemake.output.report

# Get required prefix parameter
prefix = snakemake.params.get("prefix")
if not prefix:
    raise ValueError("'prefix' parameter is required - should point to MultiQC installation directory")

# Get optional parameters
config = snakemake.params.get("config", "")
template = snakemake.params.get("template", "default")
extra_args = snakemake.params.get("extra_args", "")

# Handle logging
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Set up path to include the prefix bin directory
path_cmd = f"export PATH={prefix}/bin:$PATH && "

# Create MultiQC command
cmd = "multiqc"

# Add input files/directories
for input_file in input_files:
    cmd += f" {input_file}"

# Add output
cmd += f" -o {os.path.dirname(output_report)}"
cmd += f" -n {os.path.basename(output_report)}"

# Add template if specified
if template != "default":
    cmd += f" -t {template}"

# Add config if specified
if config:
    cmd += f" -c {config}"

# Add extra args if specified
if extra_args:
    cmd += f" {extra_args}"

# Add logging
cmd += f" {log}"

# Run commands
shell(path_cmd + cmd) 