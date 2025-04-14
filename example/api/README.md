# Biopython Translation Wrapper

A wrapper for Biopython's translation functionality, providing a simple interface for translating DNA/RNA sequences to protein sequences.

## Version

Current version: 1.1.0

For a detailed list of changes, see the [CHANGELOG.md](CHANGELOG.md).

## Usage

For general information on how to use wrappers in Snakemake, please refer to the [root README.md](../../../README.md).

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
    log:
        "proteins.log"
    wrapper:
        "file:path/to/damlab-wrappers/example/api"
```

## Input
* FASTA/FASTQ file containing DNA/RNA sequences

## Output
* FASTA file containing translated protein sequences

## Parameters

- `frame` (int, optional): Reading frame (0, 1, or 2). Defaults to 0.
- `table` (int, optional): Translation table number. Defaults to 1 (standard).
- `stop_symbol` (str, optional): Symbol used for stop codons. Defaults to "*".
- `to_stop` (bool, optional): Whether to translate only up to the first stop codon. Defaults to False.
- `cds` (bool, optional): Whether to treat the sequence as a complete CDS. Defaults to False.


## Error Handling

The wrapper includes comprehensive error handling for:
- Invalid sequences
- Invalid frame values
- Invalid translation table numbers
- Invalid stop symbols
- Translation errors

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Author
* Example Author

## Software Requirements
* [Biopython](https://biopython.org/) (tested with v1.81) 