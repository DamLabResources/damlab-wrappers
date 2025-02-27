# Generic MultiQC Log Generator

Wrapper for generating generic MultiQC-compatible placeholder log files. This wrapper is particularly useful for creating placeholder reports for samples that failed upstream QC or processing steps.

## Input
* No input files required

## Output
* YAML format log file containing all provided parameters

## Parameters
* Any parameters passed will be included in the output YAML file
* `sample_name`: Optional parameter to specify sample name (defaults to output filename stem)

## Example

```python
rule generic_multiqc_log:
    output:
        "qc/{sample}.generic.yaml"  # MultiQC-compatible log file
    params:
        sample_name = "{sample}",
        status = "FAILED",
        reason = "Low quality reads",
        qc_step = "Basecalling",
        date = "2024-03-19"
    log:
        "logs/generic_multiqc/{sample}.log"
    wrapper:
        "file:path/to/damlab-wrappers/multiqc/generic"
```

## Authors
* Will Dampier, PhD

## Notes
- This wrapper generates a standardized placeholder log file that can be parsed by MultiQC
- Useful for maintaining consistent reporting when samples fail upstream processing
- All parameters passed to the wrapper will be included in the output YAML file
- Values are converted to strings in the output YAML for consistency 