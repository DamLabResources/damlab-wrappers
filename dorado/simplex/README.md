# Dorado Simplex

Wrapper for the ONT dorado basecalling tool.
This tool performs standard basecalling on POD5 files.

## Inputs
* `pod`: Path to POD5 file or directory (required)
* `reference`: Path to reference file (optional)

## Outputs
* BAM file containing basecalled reads

## Parameters
* `model`: Model to use for basecalling (default: "sup")
* `models_directory`: Path to custom models directory (optional)
* `dorado_path`: Path to dorado executable (default: "dorado")
* `gpu`: GPU device to use (default: "all")
    * Can be "all" or specific device number

## Example

```python
rule dorado_simplex:
    input:
        pod="path/to/pod5/file.pod5",
        reference="path/to/reference.fa"  # optional
    output:
        "path/to/output.bam"
    params:
        model="sup",
        gpu="all"
    threads: 8
    wrapper:
        "file://path/to/damlab-wrappers/dorado/simplex"
```

## Authors
* Will Dampier, PhD

## Software Requirements
* [DORADO](https://github.com/nanoporetech/dorado) (tested with v0.9.1) 