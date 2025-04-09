# Biopython Translate Wrapper

This wrapper demonstrates how to create a Snakemake wrapper for a Python API. It wraps the Biopython `Seq.translate()` method to translate DNA/RNA sequences to protein sequences.

## Input
* FASTA/FASTQ file containing DNA/RNA sequences

## Output
* FASTA file containing translated protein sequences

## Parameters
* `frame` (optional, default: 0)
    Reading frame for translation (0, 1, or 2)
* `table` (optional, default: 1)
    Translation table number (1-16)
* `stop_symbol` (optional, default: "*")
    Symbol to use for stop codons
* `to_stop` (optional, default: false)
    Whether to translate to the first stop codon
* `cds` (optional, default: false)
    Whether the sequence is a complete CDS

## Example
```python
rule translate_sequences:
    input:
        "sequences.fasta"
    output:
        "proteins.fasta"
    params:
        frame=0,
        table=1,
        to_stop=True,
        stop_symbol="*"
    wrapper:
        "file:path/to/damlab-wrappers/example/api"
```

## Output Format
The output is a FASTA file containing the translated protein sequences.

## Author
* Example Author

## Software Requirements
* [Biopython](https://biopython.org/) (tested with v1.81) 