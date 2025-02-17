# HIVmetrics MultiQC Log Generator

Wrapper for generating MultiQC-compatible log files from HIV sequencing BAM files.

## Input
* BAM file containing aligned HIV sequences

## Output
* YAML format log file containing:
  - Total read counts
  - Mapped read counts
  - Supplementary alignment counts
  - Read length distribution histogram

## Parameters
* `max_size`: Maximum read length for histogram binning (default: 10,000)

## Example

```python
rule hivmetrics_multiqc_log:
    input:
        "results/{sample}.bam"
    output:
        "qc/{sample}.hivmetrics.yaml"  # MultiQC-compatible log file
    params:
        max_size = 10000  # Optional: customize max read length for histogram
    log:
        "logs/hivmetrics_multiqc/{sample}.log"
    wrapper:
        "file:path/to/damlab-wrappers/hivmetrics/multiqc"
```

## Authors
* Will Dampier, PhD

## Notes
- This wrapper generates a standardized log file that can be parsed by the DAMlab MultiQC plugin
- The log format captures key metrics from HIV sequencing data 