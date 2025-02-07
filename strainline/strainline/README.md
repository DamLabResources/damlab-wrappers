# Strainline

[Strainline](https://github.com/HaploKit/Strainline) is a tool for strain-aware long-read haplotype phasing.

## Input
- Long read sequences in FASTQ/FASTA format

## Output
One of the following:
* `haplotypes`
    FASTA file containing final haplotype sequences
* `directory`
    Directory containing all Strainline output files

## Params
* `prefix`
    Path to the Strainline installation directory
* `platform`
    Sequencing platform ("ont" or "pb", default: "ont")
* `extra_params`
    Additional parameters to pass to Strainline
* `threads`
    Number of threads to use for processing

## Example

```python
rule strainline_haplotypes:
    input:
        "reads.fastq"
    output:
        haplotypes="haplotypes.fasta"
    params:
        prefix="/path/to/strainline"
    threads: 8
    wrapper:
        "file:path/to/damlab-wrappers/strainline/strainline"


rule strainline_full:
    input:
        "reads.fastq"
    output:
        directory="strainline_output"
    params:
        prefix="/path/to/strainline",
        platform="pb",
        extra_params="--min_read_len 1000 --min_read_qual 10"
    threads: 16
    wrapper:
        "file:path/to/damlab-wrappers/strainline/strainline"
```

## Authors
* Will Dampier, PhD

## Software Requirements
* [Strainline](https://github.com/HaploKit/Strainline)
* daccord (included with Strainline)

## Notes
- The wrapper requires a properly installed Strainline environment with all dependencies
- The prefix parameter must point to the Strainline installation directory
- Either 'haplotypes' or 'directory' output must be specified, but not both 