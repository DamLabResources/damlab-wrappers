# Barcode Correction Wrapper

A wrapper for correcting barcodes in BAM files using UMI-tools clustering. This wrapper identifies and corrects barcode errors by clustering similar barcodes together based on sequence similarity.

## Version

Current version: 1.0.0 (First stable release)

For a detailed list of changes, see the [CHANGELOG.md](CHANGELOG.md).

## Installation

```bash
conda env create -f environment.yaml
```

## Usage

```python
rule correct_barcodes:
    input:
        "input.bam"
    output:
        "output.bam",
        "metrics.yaml"  # Optional metrics output
    params:
        in_tag="CR",        # SAM tag containing original barcode sequence
        out_tag="CB",       # SAM tag to store corrected barcode sequence
        barcode_length=34,  # Length of the barcode sequence
        mismatches=3,       # Maximum number of mismatches allowed when clustering barcodes
        sample_name="sample1",  # Optional sample name for metrics
        barcode_name="barcode1"  # Optional barcode name for metrics
    wrapper:
        "file:path/to/damlab-wrappers/barcode/correct"
```

## Parameters

- `in_tag` (str, required): SAM tag containing original barcode sequence
- `out_tag` (str, required): SAM tag to store corrected barcode sequence
- `barcode_length` (int, required): Length of the barcode sequence
- `mismatches` (int, optional): Maximum number of mismatches allowed when clustering barcodes. Defaults to 3.
- `sample_name` (str, optional): Sample name to include in metrics output
- `barcode_name` (str, optional): Barcode name to include in metrics output

## Input
* BAM file containing reads with barcode tags

## Output
* BAM file with corrected barcode tags
* YAML metrics file (optional) containing:
  - Unique barcode counts before and after correction
  - Total reads processed
  - Number of clusters formed
  - Number of barcodes corrected
  - Original and corrected barcode counts
  - Count histogram

## Error Handling

The wrapper includes error handling for:
- Invalid BAM files
- Missing or invalid tags
- Incorrect barcode lengths
- Clustering errors

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Authors
* Will Dampier, PhD

## Software Requirements
* [pysam](https://pysam.readthedocs.io/) (tested with v0.19.0)
* [UMI-tools](https://umi-tools.readthedocs.io/) (tested with v1.1.4) 