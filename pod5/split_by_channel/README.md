# POD5 Split by Channel

Splits POD5 files by channel using the [pod5](https://github.com/nanoporetech/pod5-file-format) tool.
This is useful for organizing nanopore data for downstream processing by channel.

## Input
- Directory containing POD5 files (`pod5_dir`). This will find all POD5 files recursively in the directory.
- List of POD5 files (`pod5_files`).

## Output
- Directory containing POD5 files split by channel

## Params
* `extra`
    Optional parameters that will be passed to pod5 subset command


## Example
```python
# Split POD5 files from a directory
rule split_by_channel:
    input:
        pod5_dir="path/to/pod5/directory"
    output:
        directory("split_by_channel")
    wrapper:
        "file:path/to/damlab-wrappers/pod5/split_by_channel"

```

## Authors
* Will Dampier, PhD

## Software Requirements
* [pod5](https://github.com/nanoporetech/pod5-file-format) (tested with v0.3.0) 