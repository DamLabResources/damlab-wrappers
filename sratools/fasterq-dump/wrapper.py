"""Wrapper for SRA Toolkit fasterq-dump command"""

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
input_path = snakemake.input[0] if snakemake.input else snakemake.params.accession
output_files = snakemake.output

# Get parameters with defaults
format = snakemake.params.get("format", "fastq")
bufsize = snakemake.params.get("bufsize", "1MB")
curcache = snakemake.params.get("curcache", "10MB")
mem = snakemake.params.get("mem", "100MB")
temp = snakemake.params.get("temp", ".")
threads = snakemake.params.get("threads", 6)
progress = snakemake.params.get("progress", True)
split_spot = snakemake.params.get("split_spot", False)
split_files = snakemake.params.get("split_files", False)
split_3 = snakemake.params.get("split_3", True)
concatenate_reads = snakemake.params.get("concatenate_reads", False)
force = snakemake.params.get("force", False)
skip_technical = snakemake.params.get("skip_technical", True)
min_read_len = snakemake.params.get("min_read_len", None)
extra = snakemake.params.get("extra", "")

# Build command arguments
args = []

# Add format argument
if format:
    args.append(f"--format {format}")

# Add buffer and cache arguments
if bufsize:
    args.append(f"--bufsize {bufsize}")
if curcache:
    args.append(f"--curcache {curcache}")
if mem:
    args.append(f"--mem {mem}")

# Add temp directory argument
if temp:
    args.append(f"--temp {temp}")

# Add threads argument
if threads:
    args.append(f"--threads {threads}")

# Add progress argument
if progress:
    args.append("--progress")

# Add split options
if split_spot:
    args.append("--split-spot")
if split_files:
    args.append("--split-files")
if split_3:
    args.append("--split-3")
if concatenate_reads:
    args.append("--concatenate-reads")

# Add force argument
if force:
    args.append("--force")

# Add skip technical argument
if skip_technical:
    args.append("--skip-technical")

# Add min read length argument
if min_read_len:
    args.append(f"--min-read-len {min_read_len}")

# Add any extra arguments
if extra:
    args.append(extra)

# Join all arguments
args_str = " ".join(args)

# Create log format string
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Execute the command
shell(
    "fasterq-dump "
    "{input_path} "
    "--outdir {snakemake.params.outdir} "
    "{args_str} "
    "{log}"
) 