# POD5 Convert FAST5 Wrapper

A wrapper for the [pod5](https://github.com/nanoporetech/pod5-file-format) tool that converts Oxford Nanopore FAST5 files to the newer POD5 format.

## Version

Current version: 1.0.0 (First stable release)

For a detailed list of changes, see the [CHANGELOG.md](CHANGELOG.md).

## Usage

For general information on how to use wrappers in Snakemake, please refer to the [root README.md](../../../README.md).

```python
rule convert_to_pod5:
    input:
        fast5_dir="path/to/fast5/directory"
    output:
        "output.pod5"
    params:
        extra="",  # Optional parameters for pod5 convert fast5
        version="1.0.0"  # Wrapper version
    threads: 4  # Optional: number of threads to use
    wrapper:
        "file://path/to/damlab-wrappers/pod5/convert_fast5"
```

## Parameters

### Input Options
- `fast5_dir` (str): Directory containing FAST5 files. The wrapper will find all FAST5 files recursively in this directory.
- `fast5_files` (str): Space-separated list of FAST5 files to convert.

### Output
- `output.pod5` (str): Path to the output POD5 file.

### Optional Parameters
- `extra` (str, optional): Additional parameters to pass to the pod5 convert fast5 command. Defaults to "".
- `version` (str, optional): Specify a version of the wrapper to use.
- `threads` (int, optional): Number of threads to use for conversion. Defaults to the number of threads specified in the Snakemake rule.

## Input
* FAST5 files (either as a directory or as a list of files)
  - Must be valid Oxford Nanopore FAST5 files
  - Can be in a single directory or distributed across subdirectories

## Output
* POD5 file
  - Contains all reads from input FAST5 files
  - Uses the newer POD5 format for improved performance

## Error Handling

The wrapper includes comprehensive error handling for:
- Missing input files or directories
- Empty directories with no FAST5 files
- Invalid output paths
- Command execution errors

## Notes
- The POD5 format is more efficient than FAST5 for downstream analysis
- Conversion can be memory-intensive for large datasets
- Using multiple threads can significantly speed up conversion

## Author
* Will Dampier, PhD

## Software Requirements
* [pod5](https://github.com/nanoporetech/pod5-file-format) (tested with v0.3.0)
* Python >= 3.10

## License
This project is licensed under the MIT License - see the LICENSE file for details. 