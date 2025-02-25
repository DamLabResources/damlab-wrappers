# Pileup Depth Calculator

A wrapper for calculating the pileup depth statistics and base composition entropy from BAM files.

## Input
* BAM file containing aligned reads (must be coordinate sorted)

## Output
* Tab-separated values (TSV) file containing:
  - Header row with column names
  - One row per covered position containing:
    - Position: Reference position
    - Depth: Number of non-gap bases
    - Entropy: Shannon entropy of base frequencies (including gaps)
    - Base counts: Count of A, C, G, T, N, and gaps (-) at this position

## Parameters
* `region` (optional)
    Region to analyze (format: "chr:start-end")
* `min_mapq` (optional, default: 0)
    Minimum mapping quality for reads to include

## Example
```python
rule calculate_pileup_depth:
    input:
        "sorted.bam"
    output:
        "pileup_stats.tsv"
    params:
        region="chr1:1000-2000",
        min_mapq=20
    wrapper:
        "file:path/to/damlab-wrappers/cigarmath/pileup"
```

## Output Format Example
```
Position    Depth    Entropy    A    C    G    T    N    -
1000        30      0.971      25   3    2    0    0    2
1001        28      0.795      0    8    0    20   0    1
1002        35      0.451      30   0    5    0    0    0
```

The entropy value ranges from 0 (all bases the same) to ~2.58 (equal distribution of all bases including gaps). Higher values indicate more base diversity at that position.

## Authors
* Will Dampier, PhD

## Software Requirements
* [pysam](https://pysam.readthedocs.io/)
* [cigarmath](https://github.com/DamLabResources/cigarmath) 