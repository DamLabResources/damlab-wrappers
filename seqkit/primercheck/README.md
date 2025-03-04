# Seqkit Primercheck

This wrapper uses seqkit amplicon to check primers against sequence reads and produces a summary CSV file.
For each read, it reports the amplicon length for each primer pair (or None if no match).

## Input
* `reads`: Input sequence file (FASTA/Q)
* `primers`: Tab-delimited primer file with columns:
    1. Primer name
    2. Forward primer (5'-3')
    3. Reverse primer (5'-3')

## Output
* Required:
    - CSV file with:
        - Rows: read IDs
        - Columns: primer names
        - Values: amplicon lengths or "None" if no match

* Optional:
    - Summary YAML file with:
        - Total number of sequences checked
        - Number of hits for each primer
        - Pairwise matrix of dual hits between primers

## Params
* `extra`: Optional parameters passed to seqkit amplicon

## Example
```python
# Basic usage
rule check_primers:
    input:
        reads="reads.fastq",
        primers="primers.txt"
    output:
        "primer_results.csv"
    wrapper:
        "file:path/to/damlab-wrappers/seqkit/primercheck"

# With summary statistics
rule check_primers_with_summary:
    input:
        reads="reads.fastq",
        primers="primers.txt"
    output:
        "primer_results.csv",
        summary="primer_summary.yaml"
    wrapper:
        "file:path/to/damlab-wrappers/seqkit/primercheck"
```

## Authors
* Will Dampier, PhD

## Software Requirements
* [seqkit](https://bioinf.shenwei.me/seqkit/) (tested with v2.3.1)
* Python packages:
    - pyyaml 