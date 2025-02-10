# Strainline MultiQC Log Generator

Wrapper for generating MultiQC-compatible log files from Strainline output directories.

## Input
* Strainline output directory containing:
  - `filter_by_abun/haplotypes.final.fa`
  - Clustering output files
  - Quality filtering logs

## Output
* YAML format log file containing:
  - Haplotype statistics (count, N50, total bases)
  - Clustering statistics
  - Quality metrics

## Example

```python
rule strainline_multiqc_log:
    input:
        "strainline_output"  # Directory containing Strainline results
    output:
        "strainline_multiqc.yaml"  # MultiQC-compatible log file
    log:
        "logs/strainline_multiqc.log"
    wrapper:
        "file:path/to/damlab-wrappers/strainline/multiqc"
```

## Authors
* Will Dampier, PhD

## Notes
- This wrapper generates a standardized log file that can be parsed by the DAMlab MultiQC plugin
- The log format is designed to capture key metrics from Strainline runs
- Future versions will expand the metrics collected based on MultiQC visualization needs 