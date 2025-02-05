"""Wrapper for MUSCLE"""

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
input_seqs = snakemake.input[0]
output_aln = snakemake.output[0]

# Get optional parameters
extra = snakemake.params.get("extra", "")

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# muscle -align input.fasta -output aligned.fasta
shell(
    "muscle -align {input_seqs} -output {output_aln} {extra} {log}"
) 