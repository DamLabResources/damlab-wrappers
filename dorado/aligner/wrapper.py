"""Wrapper for dorado aligner"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2024"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"

from snakemake.shell import shell # type: ignore
import shlex
from os import path

if "snakemake" not in locals():
    #  Keeps linters happy but doesn't impact function
    import snakemake # type: ignore

# Extract required arguments
index = snakemake.input.index
calls = snakemake.input.calls

# Separate minimap2 and dorado parameters
MM2_PARAM_FLAGS = {
    'k': '-k',
    'w': '-w',
    'I': '-I',
    'N': '-N',
    'r': '-r',
    'x': '-x',
    'junc_bed': '--junc-bed'
}

DORADO_PARAM_FLAGS = {
    'bed_file': '--bed-file',
    'max_reads': '--max-reads',
    'output_dir': '--output-dir'
}

# Build minimap2 parameter arguments
mm2_args = []
for param, flag in MM2_PARAM_FLAGS.items():
    value = snakemake.params.get(param)
    if value is not None:
        mm2_args.append(f"{flag} {shlex.quote(str(value))}")

# Handle minimap2 boolean flags
if snakemake.params.get('secondary', False):
    mm2_args.append('--secondary')
if snakemake.params.get('Y', False):
    mm2_args.append('-Y')

# Get any additional minimap2 options
extra_mm2_opts = snakemake.params.get('mm2_opts', '')
if extra_mm2_opts:
    mm2_args.append(extra_mm2_opts)

# Build dorado parameter arguments
dorado_args = []
for param, flag in DORADO_PARAM_FLAGS.items():
    value = snakemake.params.get(param)
    if value is not None:
        dorado_args.append(f"{flag} {shlex.quote(str(value))}")

# Handle dorado boolean flags
if snakemake.params.get('recursive', False):
    dorado_args.append('--recursive')
if snakemake.params.get('emit_summary', False):
    dorado_args.append('--emit-summary')
if snakemake.params.get('verbose', False):
    dorado_args.append('--verbose')

# Get dorado path
dorado_path = snakemake.params.get('dorado_path', 'dorado')

# Handle threads
threads = snakemake.threads
if threads > 0:
    dorado_args.append(f"--threads {threads}")

# Handle logging
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Combine minimap2 options if any exist
mm2_opts_str = f"--mm2-opts '{' '.join(mm2_args)}'" if mm2_args else ""

# Build and execute command
shell(
    f"{dorado_path} aligner"
    f" {' '.join(dorado_args)}"
    f" {mm2_opts_str}"
    f" {index}"
    f" {calls}"
    f" > {snakemake.output[0]}"
    f" {log}"
) 