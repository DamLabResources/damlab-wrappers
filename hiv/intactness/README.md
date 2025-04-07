# HIV Sequence Intactness Calculator

A wrapper for calculating statistics about HIV sequence intactness based on sequence lengths.

## Input
* FASTA/FASTQ file or BAM/SAM file containing HIV sequences

## Output
* YAML file containing sequence statistics:
  - total_reads: Total number of sequences analyzed
  - countable_reads: Number of sequences within countable length range
  - intact_reads: Number of sequences within intact length range
  - percent_countable: Percentage of total reads that are countable
  - percent_intact: Percentage of countable reads that are intact
  - sample_name: Name of the sample (if provided)

## Parameters
* `min_countable` (required)
    Minimum length for a sequence to be considered countable
* `max_countable` (required)
    Maximum length for a sequence to be considered countable
* `min_intact` (required)
    Minimum length for a sequence to be considered intact
* `max_intact` (required)
    Maximum length for a sequence to be considered intact
* `sample_name` (optional)
    Name of the sample to include in output

## Example
```python
rule calculate_intactness:
    input:
        "sequences.fastq"  # or .fasta, .bam, .sam
    output:
        "intactness_stats.yaml"
    params:
        min_countable=200,
        max_countable=10000,
        min_intact=8000,
        max_intact=9700,
        sample_name="patient1"
    wrapper:
        "file:path/to/damlab-wrappers/hiv/intactness"
```

## Authors
* Will Dampier, PhD

## Software Requirements
* [pysam](https://pysam.readthedocs.io/)
* [Biopython](https://biopython.org/)
* [PyYAML](https://pyyaml.org/) 