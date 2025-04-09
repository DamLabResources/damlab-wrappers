# Sequence Region Slicer

A wrapper for extracting sequences from aligned reads that overlap a specified genomic region. This wrapper uses the cigarmath library to accurately slice sequences based on CIGAR strings.

## Version

Current version: 1.0.0 (First stable release)

For a detailed list of changes, see the [CHANGELOG.md](CHANGELOG.md).

## Installation

```bash
conda env create -f environment.yaml
```

## Usage

```python
rule slice_region:
    input:
        "sorted.bam"
    output:
        fastq="region_sequences.fastq",
        metrics="region_metrics.yaml"
    params:
        region="chr1:1000-2000",
        region_name="exon1",
        sample_name="patient1",
        min_mapq=20,
        append_region_to_read_id=True
    wrapper:
        "file:path/to/damlab-wrappers/cigarmath/slice"
```

## Parameters

- `region` (str, required): Region to extract sequences from (format: "chr:start-end" or just "chr" for the entire chromosome)
- `region_name` (str, optional): A user-friendly name for the region, used in read names and metrics. Defaults to the same as region.
- `sample_name` (str, optional): Sample name to include in metrics. Defaults to an empty string.
- `min_mapq` (int, optional): Minimum mapping quality for reads to include. Defaults to 0.
- `append_region_to_read_id` (bool, optional): Whether to append region information to read IDs in the output FASTQ. Defaults to False.

## Input
* BAM file containing aligned reads (coordinate sorted recommended for performance)

## Output
* FASTQ file containing sliced sequences from the specified region
* YAML metrics file (optional) with statistics about the extraction process, including:
  - Region information
  - Sample name
  - Total segments processed
  - Segments overlapping the region
  - Whether BAM index was used
  - Whether region information was appended to read IDs

## Output Format
The output FASTQ file contains sequences extracted from the specified region, with read names formatted as:
- If `append_region_to_read_id` is True: `{original_read_name}_{region_name}_{chr}:{start}-{end}`
- If `append_region_to_read_id` is False: `{original_read_name}`

## Error Handling

The wrapper includes error handling for:
- Invalid region formats
- Missing BAM index files
- Invalid CIGAR strings
- Sequence extraction errors

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Authors
* Will Dampier, PhD

## Software Requirements
* [pysam](https://pysam.readthedocs.io/) (tested with v0.19.0)
* [cigarmath](https://github.com/DamLabResources/cigarmath) (tested with v0.1.0) 