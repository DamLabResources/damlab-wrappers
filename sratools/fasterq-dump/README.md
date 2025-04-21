# SRA Toolkit Fasterq-dump Wrapper

Version: 1.0.0

This wrapper provides a Snakemake interface for the SRA Toolkit's `fasterq-dump` command, which converts SRA files to FASTQ format.

## Input
* SRA file or directory containing SRA files

## Output
* FASTQ file(s) containing the sequence data

## Parameters
* `format` (optional, default: "fastq")
    Output format: "fastq" or "special"
* `bufsize` (optional, default: "1MB")
    Size of file buffer
* `curcache` (optional, default: "10MB")
    Size of cursor cache
* `mem` (optional, default: "100MB")
    Memory limit for sorting
* `temp` (optional, default: ".")
    Directory for temporary files
* `threads` (optional, default: 6)
    Number of threads to use
* `progress` (optional, default: true)
    Show progress
* `split_spot` (optional, default: false)
    Split spots into reads
* `split_files` (optional, default: false)
    Write reads into different files
* `split_3` (optional, default: true)
    Write single reads in special file
* `concatenate_reads` (optional, default: false)
    Write whole spots into one file
* `force` (optional, default: false)
    Force to overwrite existing file(s)
* `skip_technical` (optional, default: true)
    Skip technical reads
* `min_read_len` (optional)
    Filter by sequence length
* `extra` (optional, default: "")
    Additional parameters to pass to fasterq-dump

## Example
```python
rule convert_sra:
    input:
        "sra/SRR123456"
    output:
        "fastq/SRR123456_1.fastq",
        "fastq/SRR123456_2.fastq"
    params:
        outdir="fastq",
        split_3=True,
        threads=8,
        progress=True
    wrapper:
        "file:path/to/damlab-wrappers/sratools/fasterq-dump"
```

## Output Format
The output is one or more FASTQ files containing the sequence data. The exact number and format of output files depends on the split options used.

## Author
* William Dam

## Software Requirements
* [SRA Toolkit](https://github.com/ncbi/sra-tools) (tested with v3.0.6) 