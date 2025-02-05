"""Wrapper for FastTree"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2024"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"

from snakemake.shell import shell # type: ignore
from os import path

if "snakemake" not in locals():
    #  Keeps linters happy but doesn't impact funtion
    import snakemake # type: ignore

# Extract arguments from snakemake object
alignment = snakemake.input[0]
output = snakemake.output[0]

# Get optional parameters
extra = snakemake.params.get("extra", "")

log = snakemake.log_fmt_shell(stdout=False, stderr=True)
# FastTree -gtr -nt alignment_file > tree_file 
shell(
    "FastTree {extra} -nt {alignment} > {output} {log}"
)