# POD5 Convert FAST5

Converts Oxford Nanopore FAST5 files to the newer POD5 format using the [pod5](https://github.com/nanoporetech/pod5-file-format) tool.

## Input
- Directory containing FAST5 files (`fast5_dir`). This will find all FAST5 files recursively in the directory.
- List of FAST5 files (`fast5_files`).

## Output
- Single POD5 file

## Params
* `extra`
    Optional parameters that will be passed to pod5 convert fast5 command.

## Example
```python
rule convert_to_pod5:
    input:
        fast5_dir="path/to/fast5/directory"
    output:
        "output.pod5"
    wrapper:
        "file:path/to/damlab-wrappers/pod5/convert_fast5"

rule convert_to_pod5_with_params:
    input:
        fast5_dir="path/to/fast5/directory"
    output:
        "output.pod5"
    params:
        extra="--threads 4"
    wrapper:
        "file:path/to/damlab-wrappers/pod5/convert_fast5"
```

## Authors
* Will Dampier, PhD

## Software Requirements
* [pod5](https://github.com/nanoporetech/pod5-file-format) (tested with v0.3.0) 