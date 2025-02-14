# Dorado MultiQC Log Generator

Wrapper for generating MultiQC-compatible log files from Dorado basecalling summary files.

## Input
* `sequencing_summary.txt` file from Dorado output
  - Contains read statistics from basecalling

## Output
* YAML format log file containing:
  - Total reads
  - Mean read length
  - Mean read quality
  - Pass/fail read counts
  - Read length N50

## Example

```python
rule dorado_multiqc_log:
    input:
        "results/{sample}/sequencing_summary.txt"
    output:
        "qc/{sample}.dorado.yaml"  # MultiQC-compatible log file
    log:
        "logs/dorado_multiqc/{sample}.log"
    wrapper:
        "file:path/to/damlab-wrappers/dorado/multiqc"
```

## Authors
* Will Dampier, PhD

## Notes
- This wrapper generates a standardized log file that can be parsed by the DAMlab MultiQC plugin
- The log format captures key metrics from Dorado basecalling 