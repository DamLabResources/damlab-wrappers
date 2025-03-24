"""Wrapper for MultiQC with DAMlab plugins"""

from snakemake.shell import shell
import os

if "snakemake" not in locals():
    #  Keeps linters happy but doesn't impact funtion
    import snakemake

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
#for input_file in input_files:
#    cmd += f" {input_file}"

# Add output
basedir = os.path.dirname(output_report)
basedir = '.' if not basedir else basedir
cmd += f" -o {basedir}"
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

# Add common path of all input files
common_path = os.path.commonpath([os.path.dirname(str(f)) for f in input_files])
if common_path:
    cmd += f" {common_path}"
else:
    cmd += ' .'

# Add logging
cmd += f" {log}"

print(cmd)
# Run commands
shell(path_cmd + cmd) 