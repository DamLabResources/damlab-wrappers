# POD5 Split by Channel Wrapper

A wrapper for the [pod5](https://github.com/nanoporetech/pod5-file-format) tool that splits POD5 files by channel, organizing nanopore data for downstream processing.

## Version

Current version: 1.0.0 (First stable release)

For a detailed list of changes, see the [CHANGELOG.md](CHANGELOG.md).

## Usage

For general information on how to use wrappers in Snakemake, please refer to the [root README.md](../../../README.md).

```python
rule split_by_channel:
    input:
        pod5_dir="path/to/pod5/directory"
    output:
        directory("split_by_channel")
    params:
        extra="",  # Optional parameters for pod5 subset
        version="1.0.0"  # Wrapper version
    wrapper:
        "file://path/to/damlab-wrappers/pod5/split_by_channel"
```

## Parameters

### Input Options
- `pod5_dir` (str): Directory containing POD5 files. The wrapper will find all POD5 files recursively in this directory.
- `pod5_files` (str): Space-separated list of POD5 files to split.

### Output
- `split_by_channel` (directory): Directory containing POD5 files split by channel.

### Optional Parameters
- `extra` (str, optional): Additional parameters to pass to the pod5 subset command. Defaults to "".
- `version` (str, optional): Specify a version of the wrapper to use.

## Input
* POD5 files (either as a directory or as a list of files)
  - Must be valid POD5 files
  - Can be in a single directory or distributed across subdirectories

## Output
* Directory containing POD5 files split by channel
  - Each file contains reads from a single channel
  - Files are named according to the channel number
  - Includes a read2channel.tsv file mapping reads to channels

## Error Handling

The wrapper includes comprehensive error handling for:
- Missing input files or directories
- Empty directories with no POD5 files
- Invalid output paths
- Command execution errors

## Notes
- Splitting by channel can help organize data for parallel processing
- The read2channel.tsv file can be useful for downstream analysis
- This wrapper performs a two-step process:
  1. Creates a summary file mapping reads to channels
  2. Splits the POD5 files by channel based on this summary

## Author
* Will Dampier, PhD

## Software Requirements
* [pod5](https://github.com/nanoporetech/pod5-file-format) (tested with v0.3.0)
* Python >= 3.10

## License
This project is licensed under the MIT License - see the LICENSE file for details. 