# MUSCLE Wrapper

A wrapper for [MUSCLE](https://drive5.com/muscle/), a program for creating multiple sequence alignments of amino acid or nucleotide sequences.

## Version

Current version: 1.0.0 (First stable release)

For a detailed list of changes, see the [CHANGELOG.md](CHANGELOG.md).

## Usage

For general information on how to use wrappers in Snakemake, please refer to the [root README.md](../../../README.md).

```python
rule muscle_align:
    input:
        "sequences.fasta"
    output:
        "aligned.fasta"
    params:
        maxiters=2,
        diags=True,
        sv=False,
        refine=False,
        extra=""
    threads: 4
    wrapper:
        "file:path/to/damlab-wrappers/MSA/muscle"
```

## Parameters

- `maxiters` (int, optional): Maximum number of iterations. Defaults to None (uses MUSCLE default).
- `diags` (bool, optional): Use diagonal optimization. Defaults to False.
- `sv` (bool, optional): Use sparse vertical optimization. Defaults to False.
- `refine` (bool, optional): Refine the alignment. Defaults to False.
- `extra` (str, optional): Additional parameters to pass to MUSCLE. Defaults to "".
- `threads` (int, optional): Number of threads to use. Defaults to 1.

## Input
* FASTA file containing unaligned sequences

## Output
* FASTA file containing aligned sequences

## Error Handling

The wrapper includes comprehensive error handling for:
- Missing input files
- Invalid parameter values
- MUSCLE execution errors

## Output Format
The output is a FASTA file containing the aligned sequences, with gaps represented by "-" characters.

## Author
* Will Dampier, PhD

## Software Requirements
* [MUSCLE](https://drive5.com/muscle/) (tested with v5.1)

## License
This project is licensed under the MIT License - see the LICENSE file for details. 