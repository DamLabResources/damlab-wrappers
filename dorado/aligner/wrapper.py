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

# Get optional parameters with their flags
PARAM_FLAGS = {
    'k': '-k',
    'w': '-w',
    'I': '-I',
    'N': '-N',
    'r': '-r',
    'x': '-x',
    'junc_bed': '--junc-bed',
    'bed_file': '--bed-file',
    'max_reads': '--max-reads',
    'output_dir': '--output-dir'
}

# Build parameter arguments
param_args = []
for param, flag in PARAM_FLAGS.items():
    value = snakemake.params.get(param)
    if value is not None:
        param_args.append(f"{flag} {shlex.quote(str(value))}")

# Handle boolean flags
if snakemake.params.get('secondary', False):
    param_args.append('--secondary')
if snakemake.params.get('Y', False):
    param_args.append('-Y')
if snakemake.params.get('recursive', False):
    param_args.append('--recursive')
if snakemake.params.get('emit_summary', False):
    param_args.append('--emit-summary')
if snakemake.params.get('verbose', False):
    param_args.append('--verbose')

# Get any additional minimap2 options
mm2_opts = snakemake.params.get('mm2_opts', '')
if mm2_opts:
    param_args.append(f"--mm2-opts {shlex.quote(mm2_opts)}")

# Get dorado path
dorado_path = snakemake.params.get('dorado_path', 'dorado')

# Handle threads
threads = snakemake.threads
if threads > 0:
    param_args.append(f"--threads {threads}")

# Handle logging
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Build and execute command
shell(
    f"{dorado_path} aligner"
    f" {' '.join(param_args)}"
    f" {index}"
    f" {calls}"
    f" > {snakemake.output[0]}"
    f" {log}"
) 