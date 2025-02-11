# Strainline MultiQC Log Generator

Wrapper for generating MultiQC-compatible log files from Strainline haplotype files.

## Input
* `haplotypes.final.fa` file from Strainline output
  - Expected to be in a path like `sample_name/filter_by_abun/haplotypes.final.fa`
  - Sample name is extracted from the parent directory name

## Output
* YAML format log file containing:
  - Haplotype count
  - Haplotype length statistics (min/max/mean)
  - Haplotype frequency statistics (min/max/mean)

## Example

```python
rule strainline_multiqc_log:
    input:
        "results/{sample}/filter_by_abun/haplotypes.final.fa"
    output:
        "qc/{sample}.strainline.yaml"  # MultiQC-compatible log file
    log:
        "logs/strainline_multiqc/{sample}.log"
    wrapper:
        "file:path/to/damlab-wrappers/strainline/multiqc"
```

## Authors
* Will Dampier, PhD

## Notes
- This wrapper generates a standardized log file that can be parsed by the DAMlab MultiQC plugin
- The log format is designed to capture key metrics from Strainline haplotype reconstruction
- Sample name is extracted from the directory structure (two levels up from the haplotypes file) 