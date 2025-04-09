# Dorado Simplex

Wrapper for the ONT dorado basecalling tool.
This tool performs standard basecalling on POD5 files.

## Version

Current version: 1.0.0 (First stable release)

For a detailed list of changes, see the [CHANGELOG.md](CHANGELOG.md).

## Installation

Dorado can be installed from the Nanopore GitHub repository:

```bash
# Clone the repository
git clone https://github.com/nanoporetech/dorado.git
cd dorado

# Build dorado (requires CMake and a C++ compiler)
mkdir build && cd build
cmake ..
make -j

# Add dorado to your PATH
export PATH=$PATH:$(pwd)/bin
```

For general information on how to use wrappers in Snakemake, please refer to the [root README.md](../../../README.md).

## Usage

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

## Parameters

- `model` (str, optional): Model to use for basecalling (default: "sup")
- `models_directory` (str, optional): Path to custom models directory
- `dorado_path` (str, optional): Path to dorado executable (default: "dorado")
- `gpu` (str, optional): GPU device to use (default: "all")
  - Can be "all" or specific device number

## Input
* `pod`: Path to POD5 file or directory (required)
* `reference`: Path to reference file (optional)

## Output
* BAM file containing basecalled reads

## Error Handling

The wrapper includes error handling for:
- Missing input files
- Invalid parameters
- Dorado execution errors
- Basecalling failures

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Authors
* Will Dampier, PhD

## Software Requirements
* [DORADO](https://github.com/nanoporetech/dorado) (tested with v0.9.1) 