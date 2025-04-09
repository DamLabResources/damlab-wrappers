# Sequence Region Slicer

A wrapper for extracting sequences from aligned reads that overlap a specified genomic region.

## Input
* BAM file containing aligned reads (coordinate sorted recommended for performance)

## Output
* FASTQ file containing sliced sequences from the specified region
* Optional YAML metrics file with statistics about the extraction process

## Parameters
* `region` (required)
    Region to extract sequences from (format: "chr:start-end" or just "chr" for the entire chromosome)
* `region_name` (optional, default: same as region)
    A user-friendly name for the region, used in read names and metrics
* `sample_name` (optional, default: "")
    Sample name to include in metrics
* `min_mapq` (optional, default: 0)
    Minimum mapping quality for reads to include
* `append_region_to_read_id` (optional, default: false)
    Whether to append region information to read IDs in the output FASTQ

## Example
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

## Output Format
The output FASTQ file contains sequences extracted from the specified region, with read names formatted as:
- If `append_region_to_read_id` is True: `{original_read_name}_{region_name}_{chr}:{start}-{end}`
- If `append_region_to_read_id` is False: `{original_read_name}`

The metrics YAML file contains:
- region: The input region string
- region_name: The name of the region
- sample_name: The name of the sample
- total_segments_processed: Number of aligned segments processed from the BAM file
- segments_overlapping_region: Number of segments that overlapped the specified region
- append_region_to_read_id: Whether region information was appended to read IDs

## Authors
* Will Dampier, PhD

## Software Requirements
* [pysam](https://pysam.readthedocs.io/)
* [cigarmath](https://github.com/DamLabResources/cigarmath) 