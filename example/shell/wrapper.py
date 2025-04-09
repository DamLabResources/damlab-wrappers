"""Wrapper for seqkit translate command"""

__author__ = "Example Author"
__copyright__ = "Copyright 2024"
__email__ = "example@example.com"
__license__ = "MIT"
__version__ = "1.0.0"

from snakemake.shell import shell

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

# Get input and output files
input_file = snakemake.input[0]
output_file = snakemake.output[0]

# Get parameters with defaults
# This is a common pattern in wrappers - providing sensible defaults
frame = snakemake.params.get("frame", 0)
table = snakemake.params.get("table", 1)
trim = snakemake.params.get("trim", False)
extra = snakemake.params.get("extra", "")

# Build command arguments
# This is a common pattern - building command arguments based on parameters
args = []

# Add frame argument if specified
if frame is not None:
    args.append(f"--frame {frame}")

# Add table argument if specified
if table is not None:
    args.append(f"--table {table}")

# Add trim argument if specified
if trim:
    args.append("--trim")

# Add any extra arguments
if extra:
    args.append(extra)

# Join all arguments
args_str = " ".join(args)

# Create log format string
# This is a common pattern in wrappers - providing logging
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Execute the command
# This is the core of a shell wrapper - executing the command with the built arguments
shell(
    "seqkit translate "
    "-i {input_file} "
    "-o {output_file} "
    "{args_str} "
    "{log}"
) 