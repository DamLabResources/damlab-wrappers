# Casoff Finder

This is a tool for finding matches between gRNAs and a reference genome.

## `input_file`

The `cas-offinder` tool requires creating a custom formated configuration file that encodes:
 - A directory of fasta files
 - A search template
 - A list of gRNAs
 - The number of mismatches, insertions, and bulges.

The `casoff-finder/input_file` tool facilitates creating this file.


```python

rule:
    input:
        search_files = [],
        protospacer_file = 'list_of_protospacers.txt'
    params:
        max_mismatches = 3,
        max_dna_bulges = 2,
        max_rna_bulges = 2,
        pam = 'NGG'
    output:
        'offinder_input.txt'
    wrapper:
        ''

```

## `search`

This takes the formatted input file and an output path.
It then runs cas-offinder to search the encoded gRNAs and produces the search output.


```python

rule:
    input:
        'offinder_input.txt'
    output:
        'offinder_output.tsv'
    wrapper:
        ''
```

## `output2bed`

The `cas-offinder` output file is not easily used in downstream pipelines due to its unique format.
This wrapper uses the basic Python `csv` module to convert the output into a common bed format.


```python

rule:
    input:
        'offinder_output.tsv'
    output:
        'offinder_output.bed'
    wrapper:
        ''
```

## `multiqc`

A snakemake wrapper that creates a multiqc compatible csv file for easy summary of each gRNA.