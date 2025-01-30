# Strainline

Strainline is a tool for haplotype-aware assembly of long reads.
Due to its nature, it needs a custom environment.

The tool produces a number of output files.
You can either specify a directory, which will keep all of the files or you can specify a haplotypes output file, which will only keep the haplotype sequences. You cannot do both.


## Installation

```bash
make install PREFIX=/path/to/venv
```

This then needs to be used with the snakemake rule.

```python
rule strainline_haplotypes:
    input:
        reads = "data/reads.fastq.gz"
    output:
        haplotypes = "data/haplotypes.fasta"
    params:
        prefix = "/path/to/venv"
    wrapper:
        "file://damlab-wrappers/strainline/strainline"

rule strainline_directory:
    input:
        reads = "data/reads.fastq.gz"
    output:
        directory = directory("data/strainline_results")
    params:
        prefix = "/path/to/venv"
    wrapper:
        "file://damlab-wrappers/strainline/strainline"
```

