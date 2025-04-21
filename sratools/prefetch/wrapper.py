"""Wrapper for SRA Toolkit prefetch command"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2025"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.0.0"

from snakemake.shell import shell

if "snakemake" not in locals():
    import snakemake  # type: ignore

# Check if version is specified and compatible
if hasattr(snakemake, "params") and "version" in snakemake.params:
    requested_version = snakemake.params.version
    if requested_version != __version__:
        print(f"Warning: Requested version {requested_version} does not match wrapper version {__version__}")

# Get input and output files
accession = snakemake.params.accession
output_dir = snakemake.output[0]

# Get parameters with defaults
file_type = snakemake.params.get("type", "sra")
transport = snakemake.params.get("transport", "both")
min_size = snakemake.params.get("min_size", None)
max_size = snakemake.params.get("max_size", None)
force = snakemake.params.get("force", "no")
resume = snakemake.params.get("resume", "yes")
verify = snakemake.params.get("verify", "yes")
progress = snakemake.params.get("progress", True)
heartbeat = snakemake.params.get("heartbeat", 1)
extra = snakemake.params.get("extra", "")

# Build command arguments
args = []

# Add type argument
if file_type:
    args.append(f"--type {file_type}")

# Add transport argument
#if transport:
#    args.append(f"--transport {transport}")

# Add size limits
if min_size:
    args.append(f"--min-size {min_size}")
if max_size:
    args.append(f"--max-size {max_size}")

# Add force argument
if force:
    args.append(f"--force {force}")

# Add resume argument
if resume:
    args.append(f"--resume {resume}")

# Add verify argument
if verify:
    args.append(f"--verify {verify}")

# Add progress argument
if progress:
    args.append("--progress")

# Add heartbeat argument
if heartbeat:
    args.append(f"--heartbeat {heartbeat}")

# Add any extra arguments
if extra:
    args.append(extra)

# Join all arguments
args_str = " ".join(args)

# Create log format string
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Execute the command
shell(
    "prefetch "
    "{accession} "
    "--output-directory {output_dir} "
    "{args_str} "
    "{log}"
) 