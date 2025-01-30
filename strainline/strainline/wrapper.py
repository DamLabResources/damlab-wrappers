"""Wrapper for strainline"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2024"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"

import tempfile
import shutil
from os.path import join
import os
from snakemake.shell import shell # type: ignore


if "snakemake" not in locals():
    #  Keeps linters happy but doesn't impact funtion
    import snakemake # type: ignore

# Extract required arguments
reads = snakemake.input[0]
prefix = snakemake.params.get("prefix", "")
platform = snakemake.params.get("platform", "ont")  # default to ONT

# Determine output mode
if len(snakemake.output) != 1:
    raise ValueError("Exactly one output must be specified (either directory or haplotypes file)")

if snakemake.output.get('haplotypes', False):
    output = snakemake.output['haplotypes']
    is_directory = False
elif snakemake.output.get('directory', False):
    output = snakemake.output['directory']
    is_directory = True
else:
    raise ValueError("Exactly one output must be specified (either directory or haplotypes file)")

# Build output argument
# Get optional parameters with defaults matching help file
extra_params = snakemake.params.get("extra_params", "")

# Get number of threads
threads = f"--threads {snakemake.threads}" if snakemake.threads > 1 else ""

# Handle logging
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Ensure the strainline script path is correct
strainline_path = os.path.join(prefix, "strainline", "src", "strainline.sh")

daccord = 'daccord-0.0.10-release-20170526170720-x86_64-etch-linux-gnu'

# Build and execute command
path_cmd = f"export PATH={prefix}/bin:{prefix}/{daccord}/bin:$PATH && "
cmd = f"{strainline_path}"
cmd += f" -i {reads}"
cmd += f" -p {platform}"
cmd += f" {extra_params}"
cmd += f" {threads}"

if is_directory:
    cmd += f"-o {output}"
    shell(path_cmd + cmd)
else:
    with tempfile.TemporaryDirectory(delete=False) as tmpdir:
        cmd += f"-o {tmpdir}"
        shell(path_cmd + cmd)
        shutil.move(join(tmpdir, "filter_by_abun", "haplotypes.final.fa"), output)