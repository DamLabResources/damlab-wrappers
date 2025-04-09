# Deletion Frequency Calculator

A wrapper for calculating the frequency of reads containing deletions in specific genomic regions from BAM files. This wrapper uses the cigarmath library to accurately identify and analyze deletions in aligned reads.

## Version

Current version: 1.0.0 (First stable release)

For a detailed list of changes, see the [CHANGELOG.md](CHANGELOG.md).

## Installation

```bash
conda env create -f environment.yaml
```

## Usage

```python
rule calculate_deletion_frequency:
    input:
        "sorted.bam"
    output:
        "deletion_stats.yaml",
        "read_level_stats.csv"  # Optional CSV output
    params:
        required_region="chr1:1000-2000",
        deletion_region="chr1:1500-1600",
        region_name="exon1",
        sample_name="patient1",
        min_deletion_size=1
    wrapper:
        "file:path/to/damlab-wrappers/cigarmath/deletion_frequency"
```

## Parameters

- `required_region` (str, required): Region that reads must cover (format: "chr:start-end")
- `deletion_region` (str, required): Region where deletions are counted (format: "chr:start-end")
- `region_name` (str, optional): A user-friendly name for the region. Defaults to "region".
- `sample_name` (str, optional): Sample name to include in metrics. Defaults to "sample".
- `min_deletion_size` (int, optional): Minimum size of deletions to count. Defaults to 1.

## Input
* BAM file containing aligned reads (must be coordinate sorted)

## Output
* YAML file containing deletion frequency statistics:
  - total_reads: Total number of reads analyzed
  - reads_covering_required: Number of reads covering the required region
  - reads_with_deletion: Number of reads with deletion in specified region
  - deletion_frequency: Frequency of deletions in reads covering required region
* Optional CSV file containing read-level information:
  - read_name: Name of the read
  - covers_required: Boolean (true/false) indicating if the read covers the required region
  - has_deletion: Boolean (true/false) indicating if the read has a deletion in the specified region

## Error Handling

The wrapper includes error handling for:
- Invalid BAM files
- Invalid region formats
- Missing required parameters
- Division by zero errors

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Authors
* Will Dampier, PhD

## Software Requirements
* [pysam](https://pysam.readthedocs.io/) (tested with v0.19.0)
* [cigarmath](https://github.com/DamLabResources/cigarmath) (tested with v0.1.0)