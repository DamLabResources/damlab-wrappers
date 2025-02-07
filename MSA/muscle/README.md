# MUSCLE

[MUSCLE](https://drive5.com/muscle/) is a program for creating multiple sequence alignments of amino acid or nucleotide sequences.

## Input
- One or more sequence files in FASTA format

## Output
- Aligned sequences in FASTA format

## Params
* `extra`
    Optional parameters that will be passed to MUSCLE (e.g., "-maxiters 2")
    
## Example

```python
rule muscle_align:
    input:
        "sequences.fasta"
    output:
        "aligned.fasta"
    wrapper:
        "file:path/to/damlab-wrappers/phylo/muscle"


rule muscle_align:
    input:
        "sequences.fasta"
    output:
        "aligned.fasta"
    params:
        extra="-maxiters 2 -diags"
    wrapper:
        "file:path/to/damlab-wrappers/phylo/muscle"

```

## Authors
* Will Dampier, PhD

## Software Requirements
* [MUSCLE](https://drive5.com/muscle/) (tested with v5.1) 