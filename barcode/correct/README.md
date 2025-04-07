# Barcode Correction

A wrapper for correcting barcodes in BAM files using UMI-tools clustering.

## Input
* BAM file containing reads with barcode tags

## Output
* BAM file with corrected barcode tags

## Parameters
* `in_tag` (required)
    SAM tag containing original barcode sequence
* `out_tag` (required)
    SAM tag to store corrected barcode sequence
* `barcode_length` (required)
    Length of the barcode sequence
* `mismatches` (optional, default: 3)
    Maximum number of mismatches allowed when clustering barcodes

## Example
```python
rule correct_barcodes:
    input:
        "input.bam"
    output:
        "output.bam"
    params:
        in_tag="CR",
        out_tag="CB",
        barcode_length=34,
        mismatches=3
    wrapper:
        "file:path/to/damlab-wrappers/barcode/correct"
```

## Authors
* Will Dampier, PhD

## Software Requirements
* [pysam](https://pysam.readthedocs.io/)
* [UMI-tools](https://umi-tools.readthedocs.io/) 