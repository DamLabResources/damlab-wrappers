# SeqKit Translate Wrapper

Version: 1.0.0

This wrapper demonstrates how to create a Snakemake wrapper for a shell command. It wraps the `seqkit translate` command, which translates DNA/RNA sequences to protein sequences.

## Input
* FASTA/FASTQ file containing DNA/RNA sequences

## Output
* FASTA file containing translated protein sequences

## Parameters
* `frame` (optional, default: 0)
    Reading frame for translation (0, 1, or 2)
* `table` (optional, default: 1)
    Translation table number (1-16)
* `trim` (optional, default: false)
    Whether to trim stop codons
* `extra` (optional, default: "")
    Additional parameters to pass to seqkit translate

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
        trim=True,
        extra="--clean"
    wrapper:
        "file:path/to/damlab-wrappers/example/shell"
```

## Output Format
The output is a FASTA file containing the translated protein sequences.

## Author
* Example Author

## Software Requirements
* [seqkit](https://bioinf.shenwei.me/seqkit/) (tested with v0.16.1) 