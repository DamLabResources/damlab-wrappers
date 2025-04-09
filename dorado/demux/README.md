# Dorado Demultiplexing

Wrapper for the ONT dorado demultiplexing tool.
This tool splits reads into separate files based on the barcode sequence.

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
rule dorado_demux:
    input:
        reads="path/to/reads.bam"
    output:
        directory("path/to/output_dir")
    params:
        kit_name="SQK-RBK004",
        emit_fastq=True,
        barcode_to_output={
            "barcode01": "sample1.fastq",
            "barcode02": "sample2.fastq"
        }
    threads: 8
    wrapper:
        "file://path/to/damlab-wrappers/dorado/demux"
```

## Parameters

- `kit_name` (str, required): Name of the barcoding kit (required unless no_classify=True)
- `dorado_path` (str, optional): Path to dorado executable (default: "dorado")
- `sample_sheet` (str, optional): Path to sample sheet file
- `max_reads` (int, optional): Maximum number of reads to process
- `read_ids` (str, optional): Path to file containing read IDs to process
- `emit_fastq` (bool, optional): Output FASTQ instead of BAM (default: False)
- `emit_summary` (bool, optional): Generate summary file (default: False)
- `barcode_both_ends` (bool, optional): Require barcodes at both ends (default: False)
- `no_trim` (bool, optional): Disable barcode trimming (default: False)
- `sort_bam` (bool, optional): Sort output BAM files (default: False)
- `barcode_arrangement` (str, optional): Custom barcode arrangement
- `barcode_sequences` (str, optional): Custom barcode sequences
- `barcode_to_output` (dict, optional): Dictionary mapping barcodes to output filenames
- `tempdir` (str, optional): Directory for temporary files
- `no_classify` (bool, optional): Skip barcode classification (default: False)

## Input
* `reads`: Path to input reads file or directory (required)

## Output
* `directory`: Output directory for demultiplexed reads

## Error Handling

The wrapper includes error handling for:
- Missing input files
- Invalid parameters
- Dorado execution errors
- Demultiplexing failures

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Authors
* Will Dampier, PhD

## Software Requirements
* [DORADO](https://github.com/nanoporetech/dorado) (tested with v0.9.1)