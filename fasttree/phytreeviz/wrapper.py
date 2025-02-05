"""Wrapper for phytreeviz"""

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
input_tree = snakemake.input[0]
output_plot = snakemake.output[0]

# Get optional parameters
extra = snakemake.params.get("extra", "")
format = snakemake.params.get("format", "newick")  # Default to newick output
width = snakemake.params.get("width", 10)  # Default width in inches
height = snakemake.params.get("height", 10)  # Default height in inches

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# phytreeviz input.newick output.pdf --width 10 --height 10
shell(
    "phytreeviz -i {input_tree} -o {output_plot} --format {format} "
    "--fig_width {width} --fig_height {height} {extra} {log}"
) 