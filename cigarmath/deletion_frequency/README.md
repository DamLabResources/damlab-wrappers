# Deletion Frequency Calculator

A wrapper for calculating the frequency of reads containing deletions in specific genomic regions from BAM files.

## Input
* BAM file containing aligned reads (must be coordinate sorted)

## Output
* YAML file containing deletion frequency statistics:
  - total_reads: Total number of reads analyzed
  - reads_covering_required: Number of reads covering the required region
  - reads_with_deletion: Number of reads with deletion in specified region
  - deletion_frequency: Frequency of deletions in reads covering required region

## Parameters
* `required_region` (required)
    Region that reads must cover (format: "chr:start-end")
* `deletion_region` (required)
    Region where deletions are counted (format: "chr:start-end")

## Example
```python
rule calculate_deletion_frequency:
    input:
        "sorted.bam"
    output:
        "deletion_stats.yaml"
    params:
        required_region="chr1:1000-2000",
        deletion_region="chr1:1500-1600"
    wrapper:
        "file:path/to/damlab-wrappers/cigarmath/deletion_frequency"
```

## Authors
* Will Dampier, PhD

## Software Requirements
* [pysam](https://pysam.readthedocs.io/) 
* [cigarmath](https://github.com/DamLabResources/cigarmath)