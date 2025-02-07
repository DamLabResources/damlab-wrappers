# ClipQS

A wrapper for clipping and orienting sequences using minimap2.
This tool aligns sequences to a reference, clips unaligned regions, and ensures proper orientation.

## Input
* `sequences`
    FASTA file containing sequences to process
* `reference`
    FASTA file containing the reference sequence

## Output
- FASTA file containing processed sequences

## Params
* `min_coverage`
    Minimum coverage threshold for keeping sequences (default: 0.2)
* `include_reference`
    Whether to include the reference sequence in the output (default: True)

## Example
```python
rule clip_sequences:
    input:
        sequences="input.fasta",
        reference="reference.fasta"
    output:
        "processed.fasta"
    wrapper:
        "file:path/to/damlab-wrappers/strainline/clipqs"

rule clip_sequences:
    input:
        sequences="input.fasta",
        reference="reference.fasta"
    output:
        "processed.fasta"
    params:
        min_coverage=0.3,
        include_reference=False
    wrapper:
        "file:path/to/damlab-wrappers/strainline/clipqs"
```

## Authors
* Will Dampier, PhD

## Software Requirements
* [minimap2](https://github.com/lh3/minimap2)
* [mappy](https://pypi.org/project/mappy/) Python package 