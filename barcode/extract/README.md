# Barcode and UMI Extractor

A wrapper for extracting cellular barcodes and UMIs from BAM files using regex patterns.

## Input
* BAM file containing sequencing reads

## Output
* BAM file with added tags for barcodes and UMIs

## Parameters
* `barcode_tag` (optional, default: 'CR')
    SAM tag to use for cellular barcode sequence
* `umi_tag` (optional, default: 'OX')
    SAM tag to use for UMI sequence

Either specify a built-in pattern:
* `builtin` (optional)
    Name of built-in primer pattern set to use. Available options:
    - 'SIVmac239m2': Primers for SIVmac239m2 construct using the SIV-NFL protocol

Or specify custom patterns:
* `left_primer` (required if builtin not specified)
    Regex pattern for left primer with UMI capture group
* `right_primer` (required if builtin not specified)
    Regex pattern for right primer with UMI capture group
* `barcode_primer` (required if builtin not specified)
    Regex pattern for barcode primer with UMI capture group

## Tags Added
* CR (default): Cellular barcode sequence (34bp)
* OX (default): Combined left and right UMI sequences (36bp total)

## Example Using Built-in Pattern
```python
rule extract_barcodes:
    input:
        "input.bam"
    output:
        "output.bam"
    params:
        builtin="SIVmac239m2",
        barcode_tag="CR",
        umi_tag="OX"
    wrapper:
        "file:path/to/damlab-wrappers/barcode/extract"
```

## Example Using Custom Patterns
```python
rule extract_barcodes:
    input:
        "input.bam"
    output:
        "output.bam"
    params:
        left_primer="(PRIMER1){i<=2}(?<umi>.{18})(PRIMER2){i<=2}",
        right_primer="(PRIMER3){i<=2}(?<umi>.{18})(PRIMER4){i<=2}",
        barcode_primer="(PRIMER5){i<=2}(?<umi>.{34})(PRIMER6){i<=2}",
        barcode_tag="CR",
        umi_tag="OX"
    wrapper:
        "file:path/to/damlab-wrappers/barcode/extract"
```

## Authors
* Will Dampier, PhD

## Software Requirements
* [pysam](https://pysam.readthedocs.io/)
* [regex](https://pypi.org/project/regex/) 